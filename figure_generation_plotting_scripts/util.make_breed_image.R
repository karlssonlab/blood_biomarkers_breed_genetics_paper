library(cowplot)
library(magick)
library(tidyverse)

## NOTE: source_figures directory is not shared because it contains copyright protected images
main_dir <- "../../../../elinor_helper_code/source_figures"

# ------------------------------------------------------------------------------
# Resize breed PNGs so their heights scale with breed height (cm), preserving
# aspect ratio, and write standardized outputs to breeds_resized/.
# ------------------------------------------------------------------------------
pad_all_pngs_by_breed_height <- function() {
  src_dir  <- paste0(main_dir,"/breeds")
  out_dir  <- paste0(main_dir,"/breeds_resized")   # output folder used downstream
  min_px   <- 150
  max_px   <- 500
  
  # List all PNGs and parse breed names from filenames
  png_files <- list.files(src_dir, pattern = "\\.png$", full.names = TRUE)
  breeds_df <- tibble(
    filename = png_files,
    breed    = basename(png_files) |>
      str_remove("\\.istock\\.png$") |>
      str_replace_all("_", " ")
  )
  
  # Standardize breed naming so filenames can be joined to sheet phenotypes
  breeds_df <- standardize_breed_names(breeds_df, breed)
  
  # Join breed heights and keep only matched breeds with non-missing height
  breeds_df <- breeds_df |>
    left_join(raw_S12_BREED_PHENOS |>
                select(breed, breed.height.cm),
              by = "breed") |>
    drop_na(breed.height.cm)
  
  # Fail early if nothing matches after standardization/joining
  if (nrow(breeds_df) == 0) {
    stop("No breed filenames matched sheet breeds after standardization; check naming.")
  }
  
  # Create output directory if needed
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Map breed heights onto a target pixel height range
  breeds_df <- breeds_df |>
    mutate(scaled_height = scales::rescale(breed.height.cm, to = c(min_px, max_px)))
  
  # Compute the corresponding width for each image at the scaled height
  breeds_df <- breeds_df |>
    mutate(scaled_width = map2_dbl(filename, scaled_height, function(file, h) {
      info <- image_info(image_read(file))
      info$width * (h / info$height)
    }))
  
  # Track the maximum scaled width across breeds (useful for padding elsewhere)
  target_width <- round(max(breeds_df$scaled_width))
  
  # Resize each image to its scaled dimensions and write to breeds_resized/
  resize_and_save <- function(file, breed, new_height) {
    img  <- image_read(file)
    info <- image_info(img)
    new_width <- info$width * (new_height / info$height)
    
    # "!" forces exact dimensions (no aspect-ratio preservation by magick itself)
    img_resized <- image_resize(
      img, paste0(round(new_width), "x", round(new_height), "!")
    )
    
    # Standardized output filename convention used downstream
    out_path <- file.path(out_dir, paste0(str_replace_all(breed, " ", "_"), ".istock.png"))
    image_write(img_resized, path = out_path)
    out_path
  }
  
  # Vectorized apply across (filename, breed, scaled_height)
  pwalk(
    list(breeds_df$filename, breeds_df$breed, breeds_df$scaled_height),
    resize_and_save
  )
  
  message("Resized images saved to: ", out_dir)
}

# ------------------------------------------------------------------------------
# One-time preprocessing step: generate resized assets used by figure creation.
# ------------------------------------------------------------------------------
pad_all_pngs_by_breed_height()

# ------------------------------------------------------------------------------
# Build two one-row panels of breed images + labels for a given phenotype:
#   - low_breeds: plotset 1 (left-to-right)
#   - high_breeds: plotset 2 (left-to-right)
# Images are resized to a common max height and padded to equal canvas size,
# bottom-aligned, then combined into cowplot grids.
# ------------------------------------------------------------------------------
make_breed_figure_custom_lists <- function(
    phenotype_to_plot,
    low_breeds,
    high_breeds,
    max_px = 500,
    font_size = 5,
    img_dir = paste0(main_dir,"/breeds_resized")   # must match preprocessing output
) {
  message("Starting figure generation for phenotype: ", phenotype_to_plot)
  
  # Output directory for the padded/scaled per-phenotype images
  outdir <- file.path(paste0(main_dir,"/breeds_scaled/",phenotype_to_plot))
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Assemble requested breed list and tag which panel they belong to
  breeds <- tibble(input_breed = c(low_breeds), plotset = 1) %>%
    bind_rows(tibble(input_breed = c(high_breeds), plotset = 2)) %>%
    mutate(breed = input_breed)
  
  # Standardize breed names and build the expected PNG path in img_dir
  breeds <- standardize_breed_names(breeds, breed)
  breeds <- breeds %>%
    mutate(
      filename = file.path(img_dir, paste0(str_replace_all(breed, " ", "_"), ".istock.png"))
    )
  
  # Keep only breeds whose PNG exists in the resized image directory
  png_paths <- tibble(filename = list.files(img_dir, pattern = "\\.png$", full.names = TRUE))
  breeds <- breeds %>%
    inner_join(png_paths, by = "filename")
  
  # Fail early if none of the requested breeds have images
  if (nrow(breeds) == 0) {
    stop("No requested breeds found in '", img_dir, "'. Did you run pad_all_pngs_by_breed_height() and use the same folder?")
  }
  
  # Read original dimensions for every PNG in the directory (used to scale uniformly)
  all_breeds_info <- tibble(
    breed    = list.files(img_dir, pattern = "\\.png$") %>%
      str_replace("\\.istock\\.png$", "") %>%
      str_replace_all("_", " "),
    filename = list.files(img_dir, pattern = "\\.png$", full.names = TRUE)
  ) %>%
    mutate(
      orig_height = map_dbl(filename, ~ image_info(image_read(.x))$height),
      orig_width  = map_dbl(filename, ~ image_info(image_read(.x))$width)
    )
  
  # Use tallest image in the directory as the reference for scaling to max_px
  max_height_all <- max(all_breeds_info$orig_height)
  
  # Join dimensions for selected breeds and compute scaled sizes
  breeds <- breeds %>%
    left_join(all_breeds_info, by = c("filename", "breed")) %>%
    mutate(
      height_scaled = orig_height * (max_px / max_height_all),
      width_scaled  = orig_width  * (height_scaled / orig_height)
    )
  
  # White canvas dimensions used to pad all images to identical size
  white_width  <- max(breeds$width_scaled)
  white_height <- max_px
  
  # Output file paths mirror the input filenames but under the phenotype outdir
  breeds <- breeds %>%
   mutate(outfile = str_replace(filename, fixed(img_dir), outdir))
  
  # Resize each image then pad to a uniform white canvas, bottom-aligned ("south")
  breeds <- breeds %>%
    mutate(
      write_result = pmap(
        list(filename, width_scaled, height_scaled, outfile),
        ~ {
          img     <- image_read(..1)
          resized <- image_resize(img, geometry = paste0(round(..2), "x", round(..3), "!"))
          padded  <- image_extent(
            resized,
            geometry = paste0(round(white_width), "x", round(white_height)),
            gravity  = "south",
            color    = "white"
          )
          image_write(padded, path = ..4, format = "png")
          NULL
        }
      )
    ) %>%
    select(-write_result)
  
  # Create multi-line labels and equalize label block height by padding newlines
  breeds <- breeds %>%
    mutate(
      rel_width = 1,
      breed_label = str_replace_all(
        breed,
        "(?<!pit) (?!bull)",  # split on spaces except within the phrase "pit bull"
        "\n"
      ),
      n_lines = str_count(breed_label, "\n") + 1
    )
  max_lines <- max(breeds$n_lines)
  breeds <- breeds %>%
    mutate(breed_label = paste0(breed_label, str_dup("\n", max_lines - n_lines)))
  
  # Build a single-row cowplot grid of (image over label) grobs for a subset
  make_plot_from_breeds <- function(breeds_subset) {
    image_label_grobs <- pmap(
      list(breeds_subset$outfile, breeds_subset$breed_label),
      function(png_path, label) {
        img <- image_read(png_path)
        
        image_grob <- ggdraw() +
          draw_image(img, x = 0.5, hjust = 0.5, y = 0, vjust = 0)
        
        label_grob <- ggdraw() +
          draw_label(label, size = font_size, y = 1, vjust = 1, hjust = 0.5, lineheight = 0.8)
        
        plot_grid(image_grob, label_grob, ncol = 1, rel_heights = c(1, 0.15))
      }
    )
    
    plot_grid(
      plotlist   = image_label_grobs,
      nrow       = 1,
      align      = "v",
      axis       = "tb",
      rel_widths = breeds_subset$rel_width
    )
  }
  
  # Generate separate panels for low/high lists
  plot_low  <- make_plot_from_breeds(filter(breeds, plotset == 1))
  plot_high <- make_plot_from_breeds(filter(breeds, plotset == 2))
  
  # Return both panels for downstream composition/saving
  combined_plots <- list(low = plot_low, high = plot_high)
  combined_plots
}
