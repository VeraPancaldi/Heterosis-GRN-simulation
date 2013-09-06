#most of this code has been copied from the BoolNet function loadNetwork. 
#However, this function does not open a file to retrieve the network information,
#but uses a list that is already in the correct format, which is handed down from build.ntwk
toNet=function (tmp) 
{
    op <- c("!", "&", "\\|", "\\(", "\\)")

    targets <- sapply(tmp, function(rule) rule[1])
    factors <- sapply(tmp, function(rule) rule[2])
    probabilities <- sapply(tmp, function(rule) {
        if (length(rule) >= 3) 
            as.numeric(rule[3])
        else 1
    })
    factors.tmp <- lapply(factors, function(x) {
        sapply(op, function(y) {
            x <<- gsub(y, " ", x)
        })
        tmp <- strsplit(x, " ")[[1]]
        tmp <- unique(tmp[tmp != ""])
    })
    genes <- unique(c(targets, unname(unlist(factors.tmp))))
    isProbabilistic <- (length(unique(targets)) < length(targets))
    suppressWarnings(genes <- genes[is.na(as.integer(genes))])
    fixed <- rep(-1, length(genes))
    names(fixed) <- genes
    interactions <- list()
    for (i in 1:length(targets)) {
        target <- targets[i]
        inputGenes <- factors.tmp[[i]]
        interaction <- list()
        if (suppressWarnings(is.na(as.integer(inputGenes[1])))) {
            inputIndices <- match(inputGenes, genes)
            exp <- as.matrix(allcombn(2, length(inputIndices)) - 
                1)
            for (j in 1:length(inputIndices)) {
                tmp1 <- paste(exp[, j], collapse = ",")
                eval(parse(text = paste(inputGenes[j], "=c(", 
                  tmp1, ")", sep = "")))
            }
            tryCatch(interaction <- list(input = inputIndices, 
                func = as.numeric(eval(parse(text = factors[i]))), 
                expression = factors[i]), error = function(err) stop(paste("An error was detected in expression \"", 
                factors[i], "\": \n", err, sep = "")))
        }
        else {
            if (!isProbabilistic) 
                fixed[target] <- as.integer(inputGenes)
            interaction <- list(input = 0, func = as.integer(inputGenes), 
                expression = inputGenes)
        }
        if (isProbabilistic) {
            interaction$probability <- probabilities[i]
            interactions[[target]][[length(interactions[[target]]) + 
                1]] <- interaction
        }
        else interactions[[target]] <- interaction
    }
    onlyInputs <- setdiff(genes, targets)
    if (length(onlyInputs) > 0) {
        for (gene in onlyInputs) {
            if (isProbabilistic) 
                interactions[[gene]][[1]] = list(list(input = length(interactions) + 
                  1, func = c(0, 1), expression = gene))
            else interactions[[gene]] = list(input = length(interactions) + 
                1, func = c(0, 1), expression = gene)
        }
    }
    if (isProbabilistic) {
        wrongProb <- sapply(interactions, function(interaction) abs(1 - 
            sum(sapply(interaction, function(func) func$probability))) > 
            1e-04)
        if (any(wrongProb)) 
            stop(paste("The probabilities of gene(s) ", paste(genes[wrongProb], 
                collapse = ", "), " do not sum up to 1!", sep = ""))
    }
    res <- list(interactions = interactions, genes = genes, fixed = fixed)
    if (isProbabilistic) 
        class(res) <- c("ProbabilisticBooleanNetwork", "BooleanNetworkCollection")
    else class(res) <- "BooleanNetwork"
    return(res)
}
allcombn=function (N, n) 
{
    rownum = N^n
    sapply(n:1, function(i) {
        rep(1:N, each = N^(i - 1), len = rownum)
    })
}
