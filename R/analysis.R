#' @importFrom bio3d read.pdb
#' @export
bio3d::read.pdb


#' Analyze Model Quality
#'
#' @param pdb pdb object from read.pdb
#' @return Quality assessment results
#' @importFrom bio3d torsion.pdb
#' @export
analyze_model_quality <- function(pdb) {
  # Basic analysis
  analysis <- list(
    atoms = nrow(pdb$atom),
    residues = length(unique(pdb$atom$resno)),
    chains = unique(pdb$atom$chain),
    resolution = NA # Homology modeling usually doesn't have resolution
  )

  # Calculate Ramachandran plot statistics
  torsion <- torsion.pdb(pdb)
  if (!is.null(torsion$phi) && !is.null(torsion$psi)) {
    # Remove NA values
    phi <- torsion$phi[!is.na(torsion$phi)]
    psi <- torsion$psi[!is.na(torsion$psi)]

    if (length(phi) > 0 && length(psi) > 0) {
      # Simple Ramachandran analysis
      analysis$rama_favored <- sum(abs(phi) < 90 & abs(psi) < 90) /
        length(phi) *
        100
      analysis$rama_allowed <- sum(
        (abs(phi) >= 90 & abs(phi) < 120) | (abs(psi) >= 90 & abs(psi) < 120)
      ) /
        length(phi) *
        100
      analysis$rama_outliers <- 100 -
        analysis$rama_favored -
        analysis$rama_allowed
    }
  }

  return(analysis)
}

#' Generate Quality Report
#'
#' @param pdb_file Path to PDB file
#' @param output_dir Output directory for report
#' @return Report file path
#' @export
generate_quality_report <- function(
  pdb_file,
  output_dir = "./quality_reports"
) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  pdb <- read.pdb(pdb_file)
  quality <- analyze_model_quality(pdb)

  # Create a simple text report
  report_file <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(basename(pdb_file)), "_quality_report.txt")
  )

  report_content <- paste(
    "SWISS-MODEL Quality Report",
    "==========================",
    paste("File:", quality$file),
    paste("Atoms:", quality$atoms),
    paste("Residues:", quality$residues),
    paste("Chains:", paste(quality$chains, collapse = ", ")),
    paste("Ramachandran Favored (%):", round(quality$rama_favored, 2)),
    paste("Ramachandran Allowed (%):", round(quality$rama_allowed, 2)),
    paste("Ramachandran Outliers (%):", round(quality$rama_outliers, 2)),
    sep = "\n"
  )

  writeLines(report_content, report_file)
  message("Quality report generated: ", report_file)

  return(report_file)
}

#' Plot Ramachandran Diagram
#'
#' @param pdb pdb object from read.pdb
#' @return Plot object
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_density_2d
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ylim
#' @importFrom rlang sym
#' @importFrom tools file_path_sans_ext
#' @export
plot_ramachandran <- function(pdb) {
  torsion <- torsion.pdb(pdb)

  # Remove NA values
  phi <- torsion$phi[!is.na(torsion$phi)]
  psi <- torsion$psi[!is.na(torsion$psi)]

  if (length(phi) == 0 || length(psi) == 0) {
    warning("No torsion angles found for Ramachandran plot")
    return(NULL)
  }

  # Create data frame for plotting
  rama_data <- data.frame(phi = phi, psi = psi)

  # Create plot
  p <- ggplot(rama_data, aes(x = !!sym("phi"), y = !!sym("psi"))) +
    geom_point(alpha = 0.6, size = 1) +
    geom_density_2d(color = "red", alpha = 0.5) +
    xlim(-180, 180) +
    ylim(-180, 180) +
    labs(
      title = "Ramachandran Plot",
      x = "Phi (degrees)",
      y = "Psi (degrees)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5)
    )

  return(p)
}

#' Protein Structure Information
#'
#' @param pdb pdb object from read.pdb
#' @return Structure information
#' @export
pdb_info <- function(pdb) {
  message("Protein Structure Information:")
  message("  - Atoms: ", nrow(pdb$atom))
  message("  - Residues: ", length(unique(pdb$atom$resno)))
  message("  - Chains: ", paste(unique(pdb$atom$chain), collapse = ", "))
  # message("  - Model dimensions (Å):")
  message("  - Model dimensions (<c3><85>):")
  message("    X: ", round(max(pdb$atom$x) - min(pdb$atom$x), 2))
  message("    Y: ", round(max(pdb$atom$y) - min(pdb$atom$y), 2))
  message("    Z: ", round(max(pdb$atom$z) - min(pdb$atom$z), 2))

  return(list(
    atoms = nrow(pdb$atom),
    residues = length(unique(pdb$atom$resno)),
    chains = unique(pdb$atom$chain),
    dimensions = c(
      x = max(pdb$atom$x) - min(pdb$atom$x),
      y = max(pdb$atom$y) - min(pdb$atom$y),
      z = max(pdb$atom$z) - min(pdb$atom$z)
    )
  ))
}


#' Create Residue Composition Plot
#'
#' @param pdb pdb object from read.pdb
#' @return Plot object
#' @export
#' @importFrom ggplot2 geom_bar
#' @importFrom stats reorder
plot_residue_composition <- function(pdb) {
  # Calculate residue frequencies by type
  residue_types <- table(pdb$atom$resid)

  # Create bar plot of residue composition
  residue_df <- data.frame(
    residue = names(residue_types),
    count = as.numeric(residue_types)
  )

  p <- ggplot(
    residue_df,
    aes(x = reorder(!!sym("residue"), -!!sym("count")), y = !!sym("count"))
  ) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    labs(
      title = "Residue Composition",

      x = "Residue Type",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )

  return(p)
}

#' Visualize 3D Structure
#'
#' @param pdb pdb object from read.pdb or path to PDB file
#' @return 3D visualization object
#' @importFrom r3dmol r3dmol
#' @importFrom r3dmol m_add_model
#' @importFrom r3dmol m_set_style
#' @importFrom r3dmol m_style_cartoon
#' @importFrom r3dmol m_zoom_to
#' @importFrom bio3d write.pdb
#' @export
plot_pdb <- function(pdb) {
  if (inherits(pdb, "pdb")) {
    pdb_file <- tempfile(fileext = ".pdb")
    write.pdb(pdb, file = pdb_file)
  } else if (is.character(pdb) && file.exists(pdb)) {
    pdb_file <- pdb
  } else {
    stop("Invalid pdb input. Provide a pdb object or valid PDB file path.")
  }

  m <- r3dmol() |>
    m_add_model(pdb_file, format = 'pdb') |>
    m_set_style(style = m_style_cartoon(color = 'spectrum')) |>
    m_zoom_to()
  return(m)
}

#' Analysis Workflow
#'
#' @param pdb_file Path to PDB file
#' @param output_dir Output directory for all analysis results
#' @return Analysis results
#' @importFrom ggplot2 ggsave
#' @export
run_pdb_analysis <- function(pdb_file, output_dir = "./analysis_results") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  pdb <- read.pdb(pdb_file)

  message("Starting complete analysis workflow...")

  # 1. Quality analysis
  quality <- analyze_model_quality(pdb)

  # 2. Generate quality report
  report_file <- generate_quality_report(pdb, output_dir)

  # 3. Create Ramachandran plot
  rama_plot <- plot_ramachandran(pdb)
  rama_plot_file <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(basename(pdb_file)), "_ramachandran.png")
  )
  ggsave(rama_plot_file, rama_plot, width = 8, height = 6, dpi = 300)
  message("Ramachandran plot saved: ", rama_plot_file)

  # 4. Create residue composition plot
  residue_plot <- plot_residue_composition(pdb)
  residue_plot_file <- file.path(
    output_dir,
    paste0(
      tools::file_path_sans_ext(basename(pdb_file)),
      "_residue_composition.png"
    )
  )
  ggsave(residue_plot_file, residue_plot, width = 10, height = 6, dpi = 300)
  message("Residue composition plot saved: ", residue_plot_file)

  # 5. Text visualization
  structure_info <- pdb_info(pdb)

  results <- list(
    quality = quality,
    report_file = report_file,
    rama_plot_file = rama_plot_file,
    residue_plot_file = residue_plot_file,
    structure_info = structure_info
  )

  message("Complete analysis workflow finished!")
  return(results)
}
