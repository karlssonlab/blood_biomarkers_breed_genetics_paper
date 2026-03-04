library(tidyverse)
library(googlesheets4)
library(ggpubr)
library(cowplot)

### fix breed names 
standardize_breed_names <- function(df, column) {
  df %>%
    mutate(
      {{ column }} := tolower({{ column }}),
      {{ column }} := str_replace_all({{ column }}, "_", " "),
      {{ column }} := str_replace({{ column }}, "pitbull", "pit bull"),
      {{ column }} := str_replace({{ column }}, "st\\.", "saint"),
      {{ column }} := str_replace({{ column }}, "english bulldog", "bulldog"),
      {{ column }} := str_replace({{ column }}, "toy poodle", "poodle (toy)"),
      {{ column }} := str_replace({{ column }}, "american cocker spaniel", "cocker spaniel"),
      {{ column }} := str_replace({{ column }}, "shar pei", "chinese shar-pei"),
      {{ column }} := str_replace({{ column }}, "^entlebucher$", "entlebucher mountain dog"),
      {{ column }} := str_replace({{ column }}, "^mini ", "miniature "),
      {{ column }} := str_replace({{ column }}, "^working kelpie$", "australian kelpie"),
      {{ column }} := str_replace({{ column }}, "^kooikerhondje", "nederlandse kooikerhondje")
      
      
    )
}
### change yes/no to TRUE/FALSE
convert_yesno_to_logical <- function(df) {
  is_yesno_col <- function(x) {
    vals <- unique(tolower(na.omit(x)))
    all(vals %in% c("yes", "no", "variable"))
  }
  
  for (col_name in names(df)) {
    if (is.character(df[[col_name]]) || is.factor(df[[col_name]])) {
      if (is_yesno_col(df[[col_name]])) {
        df[[col_name]] <- sapply(tolower(as.character(df[[col_name]])), function(val) {
          if (is.na(val)) {
            NA
          } else if (val == "yes") {
            TRUE
          } else if (val == "no") {
            FALSE
          } else {
            NA
          }
        })
      }
    }
  }
  return(df)
}
library(dplyr)
library(purrr)
library(rlang)

# Helper: insert newline at separator closest to the middle (no regex)
replace_middle_break <- function(s, min_side = 4) {
  if (is.null(s) || is.na(s)) return(NA_character_)
  s <- as.character(s)
  
  # Short strings: just convert underscores to spaces
  if (nchar(s) <= 10) return(gsub("_", " ", s, fixed = TRUE))
  
  # General separator-based logic
  seps  <- c(" ", "_", "-", "\u2013", "\u2014")
  chars <- strsplit(s, "", fixed = FALSE)[[1]]
  len   <- length(chars)
  
  idx <- which(chars %in% seps)
  
  # Only allow separators that are at least `min_side` chars from each end
  valid_idx <- idx[idx > min_side & idx < (len - min_side)]
  
  # --- CASE 1: A valid natural break exists ---
  if (length(valid_idx) > 0) {
    mid <- (len + 1) / 2
    pos <- valid_idx[which.min(abs(valid_idx - mid))]
    
    sep <- chars[pos]
    keep_sep <- sep %in% c("-", "\u2013", "\u2014")
    
    left  <- if (pos > 1) paste0(chars[1:(pos - 1)], collapse = "") else ""
    right <- if (pos < len) paste0(chars[(pos + 1):len], collapse = "") else ""
    
    left  <- gsub("_", " ", left,  fixed = TRUE)
    right <- gsub("_", " ", right, fixed = TRUE)
    
    return(paste0(left, if (keep_sep) sep else "", "\n", right))
  }
  
  # --- CASE 2: No valid separator; special rule for "iso" ---
  str_len <- nchar(s)
  
  # Only consider the special "iso" rule if the string is longer than 19 characters
  if (str_len > 19) {
    # Find all "iso" positions (start index), case-insensitive
    iso_pos <- gregexpr("iso", s, ignore.case = TRUE)[[1]]
    
    if (!all(iso_pos == -1)) {
      # Break before 'iso' at position k
      # left length = k - 1, right length = str_len - (k - 1)
      # Require left > 8 and right > 8 → k > 9 and k < str_len - 7
      valid_iso <- iso_pos[iso_pos > 9 & iso_pos < (str_len - 7)]
      
      if (length(valid_iso) > 0) {
        # If multiple "iso"s, choose the one closest to the middle
        mid_str <- (str_len + 1) / 2
        pos_iso <- valid_iso[which.min(abs(valid_iso - mid_str))]
        
        left  <- substr(s, 1, pos_iso - 1)
        right <- substr(s, pos_iso, str_len)
        
        left  <- gsub("_", " ", left,  fixed = TRUE)
        right <- gsub("_", " ", right, fixed = TRUE)
        
        return(paste0(left, "-\n", right))
      }
    }
  }
  
  # --- DEFAULT: No break possible, return cleaned string ---
  return(gsub("_", " ", s, fixed = TRUE))
}

# Data frame mutator: apply to a chosen column (default: plot_label)
break_middle_label <- function(df, col = plot_label, min_side = 4) {
  col_sym <- rlang::ensym(col)
  
  df %>%
    mutate(
      {{ col_sym }} := map_chr({{ col_sym }}, ~ replace_middle_break(.x, min_side = min_side))
    )
}


format_test_result <- function(test) {
  
  if (inherits(test, "htest")) {
    
    ## Chi-squared test
    if (!is.null(test$statistic) && grepl("X-squared", names(test$statistic))) {
      stat <- unname(test$statistic)
      df   <- unname(test$parameter)
      p    <- test$p.value
      
      return(
        paste0(
          "χ² = ",
          round(stat, 2),
          ", df = ",
          df,
          ", p = ",
          format.pval(p, digits = 2, eps = 1e-300)
        )
      )
    }
    
    ## t-test
    if (!is.null(test$statistic) && grepl("t", names(test$statistic))) {
      stat <- unname(test$statistic)
      df   <- unname(test$parameter)
      p    <- test$p.value
      
      return(
        paste0(
          "t = ",
          round(stat, 2),
          ", df = ",
          round(df, 1),
          ", p = ",
          format.pval(p, digits = 2, eps = 1e-300)
        )
      )
    }
  }
  
  stop("Input must be an htest object from chisq.test() or t.test()")
}

pdf_to_plot <- function(pdf_file, density = 300, trim = TRUE,scale=1) {
  img <- image_read_pdf(pdf_file, density = density)
  if (trim) img <- image_trim(img)
  ggdraw() + draw_image(img, scale = scale)
}


