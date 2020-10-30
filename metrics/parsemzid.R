options(warn = -1)
suppressMessages(library(XML))
library(XML);
wrapper <- function() {
    context <- ''
    contextattr <- list()
    proteindata <- list(id=list(),accession=list(),length=list(),searchDatabase_ref=list())
    proteindata[["protein description"]] <- list()
    proteinattr <- NULL;
    peptidedata <- list(id=list(),sequence=list(),modifications=list())
    peptidemod <- NULL;
    themod <- NULL;
    pepevdata <- list(id=list(),peptide_ref=list(),dBSequence_ref=list(),
                      start=list(),end=list(),pre=list(),post=list(),isDecoy=list())
    pepevattr <- NULL;
    sirdata <- list(id=list(),spectrumID=list(),spectraData_ref=list())
    sirdata[["intensity of precursor ion"]] <- list()
    sirdata[["retention time"]] <- list()
    sirattr <- NULL;
    sirord <- 0;
    psmdata <- list(id=list(),rank=list(),ordinal=list(),chargeState=list(),peptide_ref=list(),experimentalMassToCharge=list(),
		    calculatedMassToCharge=list(),passThreshold=list(),spectrumIdentificationResult_ref=list())
    psmdata[["CPTAC-CDAP:PrecursorError(ppm)"]] <- list()
    psmdata[["MS-GF:QValue"]] <- list()
    psmattr <- NULL;
    counter <- 0
    start_time <- Sys.time()
    counter_tag <- "";
    DBSequence <- function(name,attrs) {
       startElement(name,attrs);
       init_count("DBSequence");
       proteinattr <<- NULL;
       for (a in names(proteindata)) {
	   if (a %in% names(attrs)) {
	       value <- attrs[[a]]
	       proteindata[[a]] <<- list(proteindata[[a]],list(value))
	       proteinattr <<- list(proteinattr,list(a))
           }
       }
       count(1000);
    };
    Peptide <- function(name,attrs) {
        startElement(name,attrs);
	init_count("Peptide");
        for (a in names(peptidedata)) {
	    if (a %in% names(attrs)) {
	       value <- attrs[[a]]
	       peptidedata[[a]] <<- list(peptidedata[[a]],list(value))
            }
	}
	peptidemod <<- list()
	count(1000)
    }
    PeptideSequence <- function(name,attrs) {
        startElement(name,attrs);
        for (a in names(peptidedata)) {
	    if (a %in% names(attrs)) {
	       value <- attrs[[a]]
	       peptidedata[[a]] <<- list(peptidedata[[a]],list(value))
            }
	}
    }
    Modification <- function(name,attrs) {
        startElement(name,attrs);
	location=as.integer(attrs[["location"]])
	monoisotopicMassDelta=as.double(attrs[["monoisotopicMassDelta"]])
	if ("residues" %in% names(attrs)) {
	    residues=attrs[["residues"]]
	} else {
	    residues="-"
	}
	themod <<- formatmod(location,monoisotopicMassDelta,residues)
    }
    formatmod <- function(location,monoisotopicMassDelta,residues) {
	s = sprintf("%d:%+.10f",location,monoisotopicMassDelta)
	s <- sub("00*$","",s)
	if (residues != '-') {
	  s = paste0(residues,s)
	}
	s
    }
    PeptideEvidence <- function(name,attrs) {
	startElement(name,attrs);
	init_count("PeptideEvidence");
	pepevattr <<- NULL;
	for (a in names(pepevdata)) {
	    if (a %in% names(attrs)) {
	       value <- attrs[[a]]
	       pepevdata[[a]] <<- list(pepevdata[[a]],list(value))
	       pepevattr <<- list(pepevattr,list(a))
            }
	}
	count(1000)
    }
    SpectrumIdentificationResult <- function(name,attrs) {
	startElement(name,attrs);
	sirattr <<- NULL;
	sirord <<- 0;
	for (a in names(sirdata)) {
	    if (a %in% names(attrs)) {
	       value <- attrs[[a]]
	       sirdata[[a]] <<- list(sirdata[[a]],list(value))
	       sirattr <<- list(sirattr,list(a))
            }
	}
    }
    SpectrumIdentificationItem <- function(name,attrs) {
	startElement(name,attrs);
	init_count("SpectrumIdentificationItem");
	psmattr <<- NULL;
	sirord <<- (sirord + 1);
	attr1 <- contextattr[["SpectrumIdentificationResult"]]
	for (a in names(psmdata)) {
	    if (a %in% names(attrs)) {
	       value <- attrs[[a]]
	       psmdata[[a]] <<- list(psmdata[[a]],list(value))
	       psmattr <<- list(psmattr,list(a))
            # } else if (a %in% names(attr1)) {
	    #    value <- attr1[[a]]
	    #    psmdata[[a]] <<- list(psmdata[[a]],list(value))
	    #    psmattr <<- list(psmattr,list(a))
            } else if (a == "spectrumIdentificationResult_ref") {
	       psmdata[[a]] <<- list(psmdata[[a]],list(attr1[["id"]]))
	       psmattr <<- list(psmattr,list(a))
	    } else if (a == "ordinal") {
	       psmdata[[a]] <<- list(psmdata[[a]],list(sirord))
	       psmattr <<- list(psmattr,list(a))
            }
	}
	count(1000)
    }
    cvParam <- function(name,attrs) {
       startElement(name,attrs);
       if (endsWith(context,".DBSequence.cvParam")) {
	   if (attrs[["name"]] %in% names(proteindata)) {
	       value <- attrs[["value"]]
	       a <- attrs[["name"]]
	       proteinattr <<- list(proteinattr,list(a))
	       proteindata[[a]] <<- list(proteindata[[a]],list(value))
           }
       } else if (endsWith(context,".Modification.cvParam")) {
	   if ("name" %in% names(attrs)) {
	       themod <<- paste0(themod,"(",attrs[["name"]],")")
	   }
       } else if (endsWith(context,".SpectrumIdentificationItem.cvParam")) {
	   if (attrs[["name"]] %in% names(psmdata)) {
	       value <- attrs[["value"]]
	       a <- attrs[["name"]]
	       psmattr <<- list(psmattr,list(a))
	       psmdata[[a]] <<- list(psmdata[[a]],list(value))
           }
       } else if (endsWith(context,".SpectrumIdentificationResult.cvParam")) {
	   if (attrs[["name"]] %in% names(sirdata)) {
	       value <- attrs[["value"]]
	       a <- attrs[["name"]]
	       sirattr <<- list(sirattr,list(a))
	       sirdata[[a]] <<- list(sirdata[[a]],list(value))
           }
       }
    };
    userParam <- function(name,attrs) {
       startElement(name,attrs);
       if (endsWith(context,".SpectrumIdentificationItem.userParam")) {
	   if (attrs[["name"]] %in% names(psmdata)) {
	       value <- attrs[["value"]]
	       a <- attrs[["name"]]
	       psmattr <<- list(psmattr,list(a))
	       psmdata[[a]] <<- list(psmdata[[a]],list(value))
           }
       } 
    };
    text <- function(thetext) {
	if (endsWith(context,".Peptide.PeptideSequence")) {
	    peptidedata[["sequence"]] <<- list(peptidedata[["sequence"]],list(thetext))
        }
    }
    startElement <- function(name,attrs) {
	context <<- paste(context,name,sep='.')
	attr <- list()
	for (n in names(attrs)) {
	    attr[[n]] = attrs[[n]]
        }
	contextattr[[name]] <<- attr
	# print(paste("BGN:",context,sep=''))
	# print(contextattr)
    };
    endElement <- function(name,attrs) {
	if (endsWith(context,".DBSequence")) {
	    proteinattr <<- unlist(proteinattr,recursive=TRUE)
	    for (a in names(proteindata)) {
		if (!(a %in% proteinattr)) {
		    proteindata[[a]] <<- list(proteindata[[a]],list(NA))
                }
            }
        } else if (endsWith(context,".Peptide")) {
	    value = paste0(unlist(peptidemod,recursive=TRUE),collapse=",")
	    peptidedata[["modifications"]] <<- list(peptidedata[["modifications"]],list(value))
        } else if (endsWith(context,".Modification")) {
	    peptidemod <<- list(peptidemod,list(themod))
	} else if (endsWith(context,".PeptideEvidence")) {
	    pepevattr <<- unlist(pepevattr,recursive=TRUE)
	    for (a in names(pepevdata)) {
		if (!(a %in% pepevattr)) {
		    pepevdata[[a]] <<- list(pepevdata[[a]],list(NA))
                }
            }
        } else if (endsWith(context,".SpectrumIdentificationItem")) {
	    psmattr <<- unlist(psmattr,recursive=TRUE)
	    for (a in names(psmdata)) {
		if (!(a %in% psmattr)) {
		    psmdata[[a]] <<- list(psmdata[[a]],list(NA))
                }
            }
	    # print(psms());
        } else if (endsWith(context,".SpectrumIdentificationResult")) {
	    sirattr <<- unlist(sirattr,recursive=TRUE)
	    for (a in names(sirdata)) {
		if (!(a %in% sirattr)) {
		    sirdata[[a]] <<- list(sirdata[[a]],list(NA))
                }
            }
        } else if (endsWith(context,".SpectrumIdentificationList")) {
	    print_count()
	}
	# print(paste("END:",context,sep=''))
	context <<- gsub(paste0("\\.",name,"$"),'',context)
	contextattr[[name]] <<- NULL;
    };
    init_count <- function(tag) {
	if (counter_tag != tag) {
	    if (counter_tag != "") {
	        print_count()
            }
	    start_time <<- Sys.time()
	    counter <<- 0
	    counter_tag <<- tag;
        }
    }
    count <- function(freq=-1) {
	counter <<- (counter + 1)
	if (freq != -1 && counter %% freq == 0) {
	    print_count();
        }
    }
    print_count <- function(final=FALSE) {
	end_time <- Sys.time()
	elapsed <- as.double(end_time - start_time,units="secs")
	# cat(sprintf("Extracted %d %s elements in %.2f seconds, %.2f per second\n",counter,counter_tag,elapsed,counter/elapsed))
	if (final) {
	    counter_tag <<- ""
        }
    }
    proteins <- function() {
	pd = list()
	for (n in names(proteindata)) {
	    pd[[n]] = unlist(proteindata[[n]],recursive=TRUE)
	    if (n == "length") {
	        pd[[n]] = as.integer(pd[[n]])
            }
	}
	data.frame(pd)
    };
    peptides <- function() {
	pd = list()
	for (n in names(peptidedata)) {
	    pd[[n]] = unlist(peptidedata[[n]],recursive=TRUE)
	}
	df1 = data.frame(pd)
    };
    peptideevidence <- function() {
	pd = list()
	for (n in names(pepevdata)) {
	    pd[[n]] = unlist(pepevdata[[n]],recursive=TRUE)
	    if (n == "start" || n == "end") {
		pd[[n]] = as.integer(pd[[n]])	
            }
	    if (n == "isDecoy") {
		pd[[n]] = sapply(pd[[n]],function(s){(s=="true")})
            }
	}
	data.frame(pd)
    };
    alignments <- function() {
	df1 = peptides()
	colnames(df1)[colnames(df1)=="id"] <- "peptide_id"

	df2 = peptideevidence()
	df2$id <- NULL

    if ("peptide_ref" %in% df2) {
	df1 = merge(df1,df2,by.x="peptide_id",by.y="peptide_ref",all.x=TRUE)

	df2 = proteins()
	colnames(df2)[colnames(df2)=="id"] <- "protein_id"

	df1 = merge(df1,df2,by.x="dBSequence_ref",by.y="protein_id",all.x=TRUE)
    }

	df1$peptide_ref <- NULL
	df1$dBSequence_ref <- NULL
	df1$protein_id <- NULL
	df1$peptide_id <- NULL
	df1$modifications <- NULL
	colnames(df1)[colnames(df1)=="sequence"] <- "peptide"
	
	unique(df1)
    };
    sii <- function() {
	pd = list()
	for (n in names(psmdata)) {
	    pd[[n]] = unlist(psmdata[[n]],recursive=TRUE)
	    if (n == "chargeState" || n == "rank" || n == "ordinal") {
		pd[[n]] = as.integer(pd[[n]])	
            }
	    if (n == "experimentalMassToCharge" || n == "calculatedMassToCharge" || n == "CPTAC-CDAP:PrecursorError(ppm)" || n == "MS-GF:QValue") {
		pd[[n]] = as.double(pd[[n]])	
            }
	    if (n == "passThreshold") {
		pd[[n]] = sapply(pd[[n]],function(s){(s=="true")})
            }
	}
	data.frame(pd)
    };
    sir <- function() {
	pd = list()
    if (length(sirdata[["id"]]) > 0) {
	for (n in names(sirdata)) {
	    pd[[n]] = unlist(sirdata[[n]],recursive=TRUE)
	    if (n == "intensity of precursor ion" || n == "retention time") {
		pd[[n]] = as.double(pd[[n]])	
            }
	    if (n == "spectrumID") {
	        pd1 <- type.convert(keyvaluetokenize(pd[[n]]))
	    }
	}
	# pure row index based merge
	# This will error if we never see spectrumID ...
	merge(data.frame(pd),pd1,by=0,all=TRUE)
    } else {
    data.frame(pd)
    }
    };
    psms <- function() {
	df1 = sii()
	df1$id <- NULL

	df2 = sir()
	colnames(df2)[colnames(df2)=="id"] <- "sir_id"

	df3 = peptides()
	colnames(df3)[colnames(df3)=="id"] <- "peptide_id"

    if ("peptide_id" %in% colnames(df3)) {
	  df1 = merge(df1,df3,by.x="peptide_ref",by.y="peptide_id",all.x=TRUE)
	  df1 = merge(df1,df2,by.x="spectrumIdentificationResult_ref",by.y="sir_id",all.x=TRUE)
    }
	df1$spectrumIdentificationResult_ref <- NULL
	df1$sir_id <- NULL
	df1$peptide_id <- NULL
	df1$peptide_ref <- NULL
	colnames(df1)[colnames(df1)=="spectraData_ref"] <- "file"
	colnames(df1)[colnames(df1)=="spectrumID"] <- "spectrum"
	df1
    };
    list(DBSequence=DBSequence,startElement=startElement,endElement=endElement,text=text,
	 cvParam=cvParam,Peptide=Peptide,PeptideSequence=PeptideSequence,Modification=Modification,
         PeptideEvidence=PeptideEvidence,SpectrumIdentificationItem=SpectrumIdentificationItem,
	 userParam=userParam,SpectrumIdentificationResult=SpectrumIdentificationResult,
	 proteins=proteins,peptides=peptides,peptideevidence=peptideevidence,psms=psms,
         sii=sii,alignments=alignments,sir=sir);
    
};
keyvaluetokenize <- function(x) {
  do.call(rbind,lapply(strsplit(x," "), function (x) { 
    x <- do.call(cbind,strsplit(x, "="))
    colnames(x) <- x[1,]
    x[-1,,drop=FALSE]}))
}
summarize <- function(name,df) {
  cat(name,"\n")
  str(df)
  summary(df)
};
parsemzid <- function(filename) {
  handlers <- wrapper();
  if (endsWith(filename,'.gz')) {
      mzidfile=gzfile(filename);
  } else {
      mzidfile=file(filename);
  }
  xmlEventParse(mzidfile,handlers=handlers);
  close(mzidfile);
  list(psms=handlers$psms(),alignments=handlers$alignments());
};

# proteins = handlers$proteins();
# summarize("Proteins",proteins)

# peptides = handlers$peptides();
# summarize("Peptides",peptides)

# peptideevidence = handlers$peptideevidence();
# summarize("PeptideEvidence",peptideevidence)

# alignments = handlers$alignments();
# summarize("Alignments",alignments)

# psms = handlers$psms();
# summarize("PSMs",psms)


