print("R script started")

rm(list=ls(all.names=TRUE))
invisible(gc())

# ==============================
# GET DIRECTORY FROM BASH
# ==============================

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
    stop("ERROR: no directory provided")
}

D <- args[1]

print(paste("Working directory:", D))

# ==============================
# LOAD FILES
# ==============================

dup_file <- file.path(D, "3D_IDduplicated.dupvar")
frq_file <- file.path(D, "3D_IDduplicated.frq")

if (file.info(dup_file)$size == 0) {

    print("No duplicated variants found")

    write.table(
        character(0),
        file.path(D, "3D_IDduplicated.exclude.txt"),
        col.names=FALSE,
        quote=FALSE,
        row.names=FALSE
    )

    quit(save="no")

}

tar_dup <- read.table(
    dup_file,
    header=TRUE,
    sep="\t",
    stringsAsFactors=FALSE
)

tar_frq <- read.table(
    frq_file,
    header=TRUE,
    stringsAsFactors=FALSE
)

# ==============================
# FUNCTION
# ==============================

fun_lowest_freq <- function(df_dup, df_frq) {

    ids_exclude <- character(0)

    if (nrow(df_dup) > 0 ) {

        for (i in 1:nrow(df_dup)) {

            ids <- unlist(strsplit(as.character(df_dup$IDS[i]), "\\s+"))

            df_ids.frq <- df_frq[df_frq$SNP %in% ids, ]

            max_id <- df_ids.frq[
                df_ids.frq$NCHROBS == max(df_ids.frq$NCHROBS),
                "SNP"
            ]

            if (length(max_id) > 1) {
                max_id <- max_id[1]
            }

            ids_exclude <- c(
                ids_exclude,
                ids[!(ids %in% max_id)]
            )
        }
    }

    return(ids_exclude)
}

# ==============================
# RUN
# ==============================

tar_exclude <- fun_lowest_freq(tar_dup, tar_frq)

# ==============================
# SAVE OUTPUT
# ==============================

out_file <- file.path(D, "3D_IDduplicated.exclude.txt")

write.table(
    tar_exclude,
    out_file,
    col.names=FALSE,
    quote=FALSE,
    row.names=FALSE
)

print(paste("Saved:", out_file))

print("R script finished")