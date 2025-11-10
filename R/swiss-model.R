# swiss-model-api.R
# Swiss-Model API Interface Functions
# Based on official documentation: https://swissmodel.expasy.org/docs/help#modelling_api

# Set SWISS-MODEL API base URL
SWISS_MODEL_BASE_URL <- "https://swissmodel.expasy.org"

#' Submit automodel project
#'
#' @param target_sequences Protein sequence string or vector of strings
#' @param project_title Project title
#' @param email Email address (optional)
#' @return Project ID
#' @export
submit_automodel <- function(
  target_sequences,
  project_title = "R Modeling Project",
  email = NULL
) {
  # Get token
  SWISS_MODEL_TOKEN <- get_swissmodel_token()
  project <- "automodel"

  # Build request body
  request_body <- list(
    target_sequences = target_sequences,
    project_title = project_title
  )

  # Add optional parameters
  if (!is.null(email)) {
    request_body$email <- email
  }

  # Send POST request
  response <- send_post_request(request_body, project = project)

  # Check response status
  check_response_status(response, project)
}


#' Submit Alignment Project
#'
#' @param target_sequences Target protein sequence
#' @param template_sequence Template protein sequence
#' @param pdb_id Template PDB ID
#' @param chain_id Chain identifier
#' @param project_title Project title
#' @param template_seqres_offset Template sequence offset (default: 0)
#' @param assembly_id Assembly ID (default: 1)
#' @return Project ID
#' @export
submit_alignment <- function(
  target_sequences,
  template_sequence,
  pdb_id,
  chain_id,
  project_title = "R Alignment Project",
  template_seqres_offset = 0,
  assembly_id = 1
) {
  # Get token
  SWISS_MODEL_TOKEN <- get_swissmodel_token()
  project <- "alignment"

  # Build request body
  request_body <- list(
    target_sequences = target_sequences,
    template_sequence = template_sequence,
    pdb_id = pdb_id,
    auth_asym_id = chain_id,
    template_seqres_offset = template_seqres_offset,
    assembly_id = assembly_id,
    project_title = project_title
  )

  # Send POST request
  response <- send_post_request(request_body, project = project)

  # Check response status
  check_response_status(response, project)
}

#' Submit User Template Project
#'
#' @param target_sequences Target protein sequence
#' @param template_coordinates_file Template coordinates file path
#' @param project_title Project title
#' @return Project ID
#' @export
submit_user_template <- function(
  target_sequences,
  template_coordinates_file,
  project_title = "R User Template Project"
) {
  # Get token
  SWISS_MODEL_TOKEN <- get_swissmodel_token()
  project <- "user_template"

  # Read template coordinates file
  if (!file.exists(template_coordinates_file)) {
    stop(
      "Template coordinates file does not exist: ",
      template_coordinates_file
    )
  }

  template_coordinates <- readLines(template_coordinates_file, warn = FALSE)
  template_coordinates <- paste(template_coordinates, collapse = "\n")

  # Build request body
  request_body <- list(
    target_sequences = target_sequences,
    template_coordinates = template_coordinates,
    project_title = project_title
  )

  # Send POST request
  response <- send_post_request(request_body, project = project)

  # Check response status
  check_response_status(response, project)
}

#' Check Project Status
#'
#' @param project_id Project ID
#' @param wait_interval Check interval in seconds (default: 10)
#' @param max_wait_time Maximum wait time in seconds (default: 1800)
#' @return Project status information
#' @export
check_project_status <- function(
  project_id,
  wait_interval = 10,
  max_wait_time = 1800
) {
  # Get token
  SWISS_MODEL_TOKEN <- get_swissmodel_token()

  start_time <- Sys.time()

  while (
    as.numeric(difftime(Sys.time(), start_time, units = "secs")) < max_wait_time
  ) {
    # Get project status
    response <- request(paste0(
      SWISS_MODEL_BASE_URL,
      "/project/",
      project_id,
      "/models/summary/"
    )) |>
      req_headers("Authorization" = paste("Token", SWISS_MODEL_TOKEN)) |>
      req_perform()

    status <- resp_status(response)
    if (status == 200) {
      status_info <- resp_body_json(response)
      status_text <- status_info$status

      elapsed_time <- as.numeric(difftime(
        Sys.time(),
        start_time,
        units = "secs"
      ))
      message(
        "Project status: ",
        status_text,
        " [Waited: ",
        round(elapsed_time),
        " seconds]"
      )

      if (status_text %in% c("COMPLETED", "FAILED")) {
        return(status_info)
      }
    } else {
      error_msg <- resp_body_string(response)
      message("Failed to get status: ", status, " - ", error_msg)
    }

    # Wait for specified interval
    Sys.sleep(wait_interval)
  }

  stop("Timeout reached (", max_wait_time, " seconds)")
}

#' Get Model Coordinates URLs
#'
#' @param status_info Project status information object
#' @param file_format File format ("pdb" or "cif", default: "pdb")
#' @return List of model coordinate URLs
#' @export
get_model_coordinates_urls <- function(status_info, file_format = "pdb") {
  if (status_info$status != "COMPLETED") {
    stop("Project not completed, current status: ", status_info$status)
  }

  model_list <- status_info$models
  url_list <- list()

  for (i in seq_along(model_list)) {
    model <- model_list[[i]]

    if (file_format == "pdb" && !is.null(model$coordinates_url)) {
      url_list[[paste0("model_", i)]] <- model$coordinates_url
    } else if (file_format == "cif" && !is.null(model$modelcif_url)) {
      url_list[[paste0("model_", i)]] <- model$modelcif_url
    }
  }

  return(url_list)
}

#' Download Model Files
#'
#' @param coordinates_url Coordinates file URL
#' @param output_file Output file path
#' @return Downloaded file path
#' @importFrom httr2 request
#' @importFrom httr2 req_headers
#' @importFrom httr2 req_perform
#' @importFrom httr2 resp_status
#' @importFrom httr2 resp_body_raw
#' @importFrom httr2 resp_body_string
#' @export
download_model_file <- function(coordinates_url, output_file) {
  # Get token
  SWISS_MODEL_TOKEN <- get_swissmodel_token()

  response <- request(coordinates_url) |>
    req_headers("Authorization" = paste("Token", SWISS_MODEL_TOKEN)) |>
    req_perform()

  status <- resp_status(response)
  if (status == 200) {
    writeBin(resp_body_raw(response), output_file)
    message("File downloaded successfully: ", output_file)
    return(output_file)
  } else {
    error_msg <- resp_body_string(response)
    stop("Download failed: ", status, " - ", error_msg)
  }
}

#' Bulk Download Project Results
#'
#' @param from_datetime Start time (optional, format: YYYY-MM-DD HH:MM:SS)
#' @param to_datetime End time (optional, format: YYYY-MM-DD HH:MM:SS)
#' @param coordinates_type Coordinates file type ("cif" or "pdb", default: "pdb")
#' @param output_dir Output directory
#' @return Downloaded file path
#' @importFrom httr2 request
#' @importFrom httr2 req_url_query
#' @importFrom httr2 req_headers
#' @importFrom httr2 req_perform
#' @importFrom httr2 resp_status
#' @importFrom httr2 resp_body_json
#' @importFrom httr2 resp_body_raw
#' @export
bulk_download_projects <- function(
  from_datetime = NULL,
  to_datetime = NULL,
  coordinates_type = "pdb",
  output_dir = "./bulk_download"
) {
  # Get token
  SWISS_MODEL_TOKEN <- get_swissmodel_token()

  # Build request parameters
  params <- list()
  if (!is.null(from_datetime)) {
    params$from_datetime <- from_datetime
  }
  if (!is.null(to_datetime)) {
    params$to_datetime <- to_datetime
  }
  if (!is.null(coordinates_type)) {
    params$coordinates_type <- coordinates_type
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Start bulk download job
  response <- request(paste0(SWISS_MODEL_BASE_URL, "/projects/download/")) |>
    req_url_query(!!!params) |>
    req_headers("Authorization" = paste("Token", SWISS_MODEL_TOKEN)) |>
    req_perform()

  status <- resp_status(response)
  if (status != 200) {
    stop("Bulk download request failed: ", status)
  }

  download_id <- resp_body_json(response)$download_id

  # Wait for bulk download to complete
  message("Waiting for bulk download to complete...")
  while (TRUE) {
    Sys.sleep(5)

    # Check status
    status_response <- request(paste0(
      SWISS_MODEL_BASE_URL,
      "/projects/download/",
      download_id,
      "/"
    )) |>
      req_headers("Authorization" = paste("Token", SWISS_MODEL_TOKEN)) |>
      req_perform()

    status_info <- resp_body_json(status_response)

    if (status_info$status %in% c('COMPLETED', 'FAILED')) {
      break
    }
  }

  if (status_info$status == 'COMPLETED') {
    # Download result file
    download_url <- status_info$download_url
    output_file <- file.path(
      output_dir,
      paste0("bulk_download_", download_id, ".zip")
    )

    download_response <- request(download_url) |> req_perform()
    status <- resp_status(download_response)
    if (status == 200) {
      writeBin(resp_body_raw(download_response), output_file)
      message("Bulk download completed: ", output_file)
      return(output_file)
    } else {
      stop("Failed to download bulk results")
    }
  } else {
    stop("Bulk download failed")
  }
}

#' Complete Automodel Workflow
#'
#' @param protein_sequence Protein sequence
#' @param project_title Project title
#' @param output_dir Output directory
#' @param wait_interval Check interval
#' @param max_wait_time Maximum wait time
#' @return Modeling results
#' @importFrom jsonlite write_json
#' @export
run_automodel_workflow <- function(
  protein_sequence,
  project_title = "R Automodel",
  output_dir = "./modeling_results",
  wait_interval = 30,
  max_wait_time = 1800
) {
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  message("Starting automodel workflow...")

  # 1. Submit modeling task
  project_id <- submit_automodel(protein_sequence, project_title)
  message("Project ID: ", project_id)

  # 2. Wait for task completion
  status_info <- check_project_status(project_id, wait_interval, max_wait_time)

  if (status_info$status == "COMPLETED") {
    message("Modeling completed successfully!")

    # 3. Get model coordinate URLs
    coordinates_urls <- get_model_coordinates_urls(status_info, "pdb")

    # 4. Download all model files
    downloaded_files <- list()
    for (model_name in names(coordinates_urls)) {
      output_file <- file.path(
        output_dir,
        paste0(project_id, "_", model_name, ".pdb")
      )
      downloaded_file <- download_model_file(
        coordinates_urls[[model_name]],
        output_file
      )
      downloaded_files[[model_name]] <- downloaded_file
    }

    # Save project information
    project_info_file <- file.path(output_dir, paste0(project_id, "_info.json"))
    write_json(status_info, project_info_file, pretty = TRUE, auto_unbox = TRUE)

    return(list(
      project_id = project_id,
      status = "COMPLETED",
      downloaded_files = downloaded_files,
      project_info_file = project_info_file
    ))
  } else {
    stop("Modeling failed, status: ", status_info$status)
  }
}
