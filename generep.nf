#!/usr/bin/env nextflow

params.resultsdir = workflow.launchDir + "/results/"

shufnames = Channel.of(1..params.nshuf)

process InitializeGenerep {

	output:
	file("real.csv") into realcleaned
	file("real.csv") into realboot

	"""
	#!/bin/bash
	mkdir -p ${params.resultsdir}
	genereptools.py repeat ${params.expressions} real.prez.csv
	genereptools.py zscore real.prez.csv real.csv
	"""
}

Channel
// Generates the names of all shuffled directories: shuffled1, shuffled2, etc
	.from shufnames
	.map { "shuffled" + it }
	.set { shufdirs }

process Randomize {

	input:
	val real from realcleaned
	val shuf from shufdirs

	output:
	tuple val(shuf), file("${shuf}.csv") into shuffled

	"""
	#!/bin/bash
	apple.py random -o ${shuf}.csv -gc ${params.genecols} ${real}
	"""
}

process BootstrapReal {

	time "2:00:00"
	memory "2G"

	input:
	val real from realboot

	output:
	file("real-*.csv") into realbootstraps

	"""
	#!/bin/bash
	cp ${real} real.csv
	apple.py bootstrap -gc ${params.genecols} real.csv ${params.nrounds}
	"""
}

process Bootstrap {
	
	time "2:00:00"
	memory "2G"

	input:
	tuple val(shuf), file(shuffile) from shuffled
	
	output:
	file("shuffled?-*.csv") into bootstraps

	"""
	#!/bin/bash
	cp ${shuffile} ${shuf}.csv
	apple.py bootstrap -gc ${params.genecols} ${shuf}.csv ${params.nrounds}
	"""
}

Channel
	.from realbootstraps
	.concat(bootstraps)
	.flatten()
	.map { file -> tuple(file.baseName.tokenize("-")[0], file.baseName, file) }
	.set { allbootstraps }

process Aracne {

	time "80:00:00"
	memory "30G"

	input:
	tuple val(key), val(basename), file(boot) from allbootstraps

	output:
	tuple val(key), file("${basename}.adj") into alladjs

	"""
	#!/bin/bash
	# Running aracne on $basename $boot  ${basename}.adj

	# Figure out where aracne2 is
	AR=\$(which aracne2)
	if [ "\$?" == "0" ];
	then
	  ARDIR=\$(dirname \$AR)
	  \$ARDIR/aracne2 -i ${boot} -o ${basename}.adj -l ${params.tfs} \
          -a adaptive_partitioning \
          -H \$ARDIR \
	  -e ${params.dpi}
	else
	  echo "Error: aracne2 not found in PATH."
	  exit 1
	fi
	"""
}

Channel
	.from alladjs
	.groupTuple()
	.into { adjgroups; adjgroups2 }

process Consensus {
	
	time "20:00:00"
	memory "20G"

	input:
	tuple val(key), file(adjs) from adjgroups

	output:
	tuple val(key), file("${key}.adj"), file("${key}.counts.csv") into origConsensus
	file("${key}.counts.csv") into allconsensus

	"""
	#!/bin/bash
	apple.py consensus -s 1 -c ${key}.counts.csv ${key}.adj ${adjs}
	"""
}

Channel
	// Take the adj files created by Consensus and collect them in a list, for BestSupport
	.from allconsensus
	.collect()
	.set { findbestsupport }

process BestSupport {
	
	executor "local"

	input:
	val(files) from findbestsupport

	output:
	file("optimal-support.csv") into bestsupport
	file("optimal-support.csv") into thresholds1

	"""
#!/bin/bash
genereptools.py suppfilter ${params.nrounds} "${files}"
	"""


//	"""
// #!/usr/bin/env python
// import genereputils
// F = genereputils.GRFilterSupport(None)
// F.runNF("${files}", ${params.nrounds})
//	"""
}

process FilterBySupport {
	// Create a new consensus filtering edges by support.

	time "10:00:00"
	memory "20G"

	input:
	tuple val(key), file(adjs) from adjgroups2
	file(optsupport) from bestsupport
	
	output:
	tuple val(key), file("${key}.supp.adj") into suppadj
	file("${key}.mi.hist.csv") into mihistogram

	"""
	#!/bin/bash
	supp=`grep OptSupport ${optsupport} | cut -f 2`
	apple.py consensus -s \$supp -c ${key}.supp.counts.csv ${key}.supp.adj ${adjs}
	apple.py histogram -n 1000 -o ${key}.mi.hist.csv ${key}.supp.adj
	cp ${key}.supp.adj ${params.resultsdir}
	"""
}

Channel
// Collect all *.mi.hist.csv files to determine optimal MI	
	.from mihistogram
	.collect()
	.set { allmihist }

process OptimalMI {
	// Analyze MI histograms to find the optimal MI threshold
	executor "local"

	input:
	val allmihist from allmihist
	
	output:
	file("optimal-mi.csv") into bestmi
	file("optimal-mi.csv") into thresholds2

	"""
#!/bin/bash
genereptools.py mifilter "${allmihist}"
"""

//	"""
// #!/usr/bin/env python
// import genereputils
// F = genereputils.GRFilterMI(None)
// F.runNF("${allmihist}")
// """
}

process FilterByMI {
	// Filter the consensus adj files using the optimal MI threshold

	time "5:00:00"
	memory "20G"

	input:
	tuple val(key), file(adj) from suppadj
	file(optmi) from bestmi
	
	output:
	tuple val(key), file("${key}.mi.adj") into filterBySumMI
	file("${key}.summi.hist.csv") into summihistogram

	"""
	#!/bin/bash
	mi=`grep OptMI ${optmi} | cut -f 2`
	apple.py filter -o ${key}.mi.adj ${adj} \$mi
	apple.py histogram -n 1000 -s -o ${key}.summi.hist.csv ${key}.mi.adj
	cp ${key}.mi.adj ${params.resultsdir}
	"""
}

Channel
// Collect all *.summi.hist.csv files to determine optimal SumMI	
	.from summihistogram
	.collect()
	.set { allsummihist }

process OptimalSumMI {
	// Analyze MI histograms to find the optimal MI threshold
	executor "local"

	input:
	val allsummihist from allsummihist
	
	output:
	file("optimal-summi.csv") into bestsummi
	file("optimal-summi.csv") into thresholds3

	"""
#!/bin/bash
genereptools.py summifilter "${allsummihist}"
	"""

// #!/usr/bin/env python
// import genereputils
// F = genereputils.GRFilterSumMI(None)
// F.runNF("${allsummihist}")

}

process FilterBySumMI {
	// Filter the consensus adj files using the optimal MI threshold

	time "5:00:00"
	memory "20G"

	input:
	tuple val(key), file(adj) from filterBySumMI
	file(optsummi) from bestsummi
	
	output:
	tuple val(key), file("${key}.summi.adj"), file("${key}.final.hist.csv") into filterFinal
	file("${key}.summi.adj") into adjstats3

	"""
	#!/bin/bash
	summi=`grep OptSumMI ${optsummi} | cut -f 2`
	apple.py filter -t -o ${key}.summi.adj ${adj} \$summi
	apple.py histogram -n 1000 -o ${key}.final.hist.csv ${key}.summi.adj
	cp ${key}.summi.adj ${params.resultsdir}
	"""
}

Channel
	.from filterFinal
	.collect()
	.set { finaladjs }

process FilterFinal {

	input:
	val(allfinaladjs) from finaladjs

	output:
	file("${params.label}.conn.csv") into connections

	"""
	#!/bin/bash

	res=${params.resultsdir}
	apple.py stats \$res/*.supp.adj \$res/*.mi.adj \$res/*.summi.adj > \$res/adj-stats.csv

	listToTable.py "${allfinaladjs}" | sort > ADJS
	realadj=`grep ^real ADJS | cut -f 2`
	realhist=`grep ^real ADJS | cut -f 3`
	shufadjs=`grep -v ^real ADJS | cut -f 2`
	shufhists=`grep -v ^real ADJS | cut -f 3`

	# Generate final histograms
	finalhist=real.vs.shuffled.avg.hist.csv
	histfiles=""
	while read LABEL ADJFILE HISTFILE;
	do 
	  if [ "\$LABEL" != "real" ];
	  then
	    outfile="real.vs.\${LABEL}.final.hist.csv"
	    histfiles="\$outfile \$histfiles"
	    genereptools.py histogram compare \$outfile \$realhist \$HISTFILE
	  fi
	done < ADJS
	genereptools.py histogram average \$finalhist \$histfiles
	
	# Convert final adj file to cytoscape
	cyto=${params.label}.cyto.csv

	rm -f unsorted.cyto
	if [ -f \$finalhist ];
	then
	  apple.py convert ac -f \$finalhist \$realadj unsorted.cyto
	else
	  apple.py convert ac \$realadj unsorted.cyto
	fi
	head -1 unsorted.cyto > \$cyto
	tail -n +2 unsorted.cyto | sort >> \$cyto

	# If translation file is provided, translate gene identifiers
	if [ -f "${params.ids}" ];
	then
	  cytopre=\$cyto
	  cyto=${params.label}.tr.csv
	  apple.py translate ${params.ids} \$cytopre \$cyto
	  adjpre=\$realadj
	  realadj=real.tr.adj
	  apple.py convert ca \$cyto \$realadj
	fi

	# Produce connections file and print top genes
	connections=${params.label}.conn.csv
	apple.py convert co \$cyto \$connections
	cut -f 1,2 \$connections | head -50 > tophubs.txt

	# Produce CX file
	cx=${params.label}.cx
	cxtop=${params.label}-tophubs.cx
	echo "${params.cxattributes}" > attributes.txt
	apple.py convert cx -a attributes.txt \$cyto \$cx
	apple.py extract -a -c -o tophubnet.cy \$cyto tophubs.txt
	apple.py convert cx -a attributes.txt tophubnet.cy \$cxtop

	# Copy output files to resultsdir
	cp \$cyto \$connections \$realadj \$finalhist \$cx \$cxtop tophubs.txt ${params.resultsdir}
	"""
}

process CollectThresholds {

	input:
	path(thr1) from thresholds1
	path(thr2) from thresholds2
	path(thr3) from thresholds3

	"""
	#!/bin/bash
	cat ${thr1} ${thr2} ${thr3} > ${params.resultsdir}/thresholds.csv
	"""
}


