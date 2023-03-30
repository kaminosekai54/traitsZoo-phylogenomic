settings ={
"path" : {
"rawData": "rawData/",
"data": "data/"
},
"referenceProteom":"Amphibalanus-amphitrite",
# "referenceProteom":"Oithona-nana",
"pathToTrimal": "~/trimAl/source/"
}

settings["path"]["renamedFasta"] = settings["path"]["data"] + "fastaFiles/"
settings["path"]["tblFiles"] = settings["path"]["data"] + "tblFiles/"
settings["path"]["indexedFiles"] =  settings["path"]["data"] + "indexedFiles/"
settings["path"]["tmpBash"] = settings["path"]["data"] + "tmpBash/"
settings["path"]["diamondMatchs"] = settings["path"]["data"] + "diamondMatchs/"
settings["path"]["RBH"] = settings["path"]["data"] + "RBH/"
settings["path"]["rawRBH"] = settings["path"]["RBH"] + "raw/"
settings["path"]["filteredRBH"] = settings["path"]["RBH"] + "filtered/"
settings["path"]["alignments"]= settings["path"]["data"] + "alignments/"
settings["path"]["rawAlignments"]= settings["path"]["alignments"] + "raw/"
settings["path"]["trimmedAlignments"]= settings["path"]["alignments"] + "trimmed/"
settings["path"]["trees"] =settings["path"]["data"] + "trees/",


def getSettings(): return settings
