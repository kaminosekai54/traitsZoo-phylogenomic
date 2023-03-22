settings ={
"path" : {
"rawData": "rawData/",
"data": "data/"
},
"referenceProteom":"Oithona-nana"
}

settings["path"]["renamedFasta"] = settings["path"]["data"] + "fastaFiles/"
settings["path"]["tblFiles"] = settings["path"]["data"] + "tblFiles/"
settings["path"]["indexedFiles"] =  settings["path"]["data"] + "indexedFiles/"
settings["path"]["tmpBash"] = settings["path"]["data"] + "tmpBash/"
settings["path"]["RBH"] = settings["path"]["data"] + "RBH/"
settings["path"]["alignments"]= settings["path"]["data"] + "alignments/"
settings["path"]["trees"] =settings["path"]["data"] + "trees/",


def getSettings(): return settings
