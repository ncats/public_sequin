# create error log file
makeErrorLog = function() {
    if (file.exists("error_log.tx")) {
        file.remove("error_log.tx")
    }
    file.create("error_log.tx")
}
makeErrorLog()

# Function to write errors to the log file
write_error_to_log <- function(message, block_name) {
    error_message <- paste(Sys.time(), "Error in", block_name, ":", message, "\n", sep = " ")
    cat(error_message, file = "error_log.tx", append = TRUE)
}

