
.onLoad <- function(libname, pkgname, where) {

    if(.Platform$OS.type == "windows" && require(Biobase) && interactive()
        && .Platform$GUI ==  "Rgui"){
        addVigs2WinMenu("Resourcerer")
    }
    require("AnnotationDbi", quietly = TRUE) ||
    stop("Package AnnotationDbi unavailable!")
}
