#!/usr/bin/env Rscript
# Read-only R/TAM runtime identity, including installed code and library bytes.
for (p in c("TAM","numDeriv","jsonlite","digest")) suppressPackageStartupMessages(loadNamespace(p))
stopifnot(as.character(packageVersion("TAM")) == "4.3.25")
packages <- sort(loadedNamespaces())
paths <- setNames(lapply(packages, function(p) system.file(package=p)), packages)
stopifnot(all(nzchar(unlist(paths))))
files <- unique(unlist(lapply(paths,function(p) list.files(p,recursive=TRUE,full.names=TRUE))))
files <- sort(files[!dir.exists(files)])
libR <- file.path(R.home("lib"),"libR.dylib")
if (file.exists(libR)) files <- sort(unique(c(files,libR)))
result <- list(R=R.version.string,platform=R.version$platform,library_paths=.libPaths(),
    packages=setNames(lapply(packages,function(p) as.character(packageVersion(p))),packages),
    file_sha256=setNames(lapply(files,function(p) digest::digest(file=p,algo="sha256")),files))
cat(jsonlite::toJSON(result,auto_unbox=TRUE,digits=NA,null="null"))
