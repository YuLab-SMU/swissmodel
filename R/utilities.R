#' Set API Token
#'
#' @param token API token string
#' @return No return value
#' @export
set_swissmodel_token <- function(token) {
  if (is.null(token) || nchar(token) == 0) {
    stop(
      "Please provide a valid API token. You can get it from your SWISS-MODEL account page: https://swissmodel.expasy.org/account"
    )
  }
  options(SWISS_MODEL_TOKEN = token)
  message("API token set successfully")
}

get_swissmodel_token <- function() {
  token <- getOption("SWISS_MODEL_TOKEN")
  if (is.null(token)) {
    stop("API token not set. Please use set_swissmodel_token() to set it.")
  }
  return(token)
}


# Send POST Request to SWISS-MODEL API
#' @importFrom httr2 request
#' @importFrom httr2 req_headers
#' @importFrom httr2 req_body_json
#' @importFrom httr2 req_perform
send_post_request <- function(request_body, project) {
  project <- match.arg(
    tolower(project),
    choices = c("automodel", "alignment", "user_template")
  )

  # Get token
  SWISS_MODEL_TOKEN <- get_swissmodel_token()

  response <- sprintf("%s/%s", SWISS_MODEL_BASE_URL, project) |>
    request() |>
    req_headers(
      "Authorization" = paste("Token", SWISS_MODEL_TOKEN),
      "Content-Type" = "application/json"
    ) |>
    req_body_json(request_body) |>
    req_perform()

  return(response)
}


# Check Response Status
#' @importFrom httr2 resp_status
#' @importFrom httr2 resp_body_json
#' @importFrom httr2 resp_body_string
check_response_status <- function(response, project) {
  status <- resp_status(response)
  if (status == 202) {
    message(
      toupper(project),
      " task submitted successfully!"
    )
    result <- resp_body_json(response)
    return(result$project_id)
  } else if (status == 200) {
    message("Project already completed/failed (same input and SMTL version)")
    result <- resp_body_json(response)
    return(result$project_id)
  } else if (status == 429) {
    stop(
      "Request rate too high. Current limits: rapid submission rate 100/minute, prolonged rate 2000/6 hours"
    )
  } else {
    error_msg <- resp_body_string(response)
    stop("Submission failed: ", status, " - ", error_msg)
  }
}
