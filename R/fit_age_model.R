
#' Fit a model (or age-length key) that will be used to estimate ages
#'
#' This function will create the appropriate model (or age-length key) that can
#' then be used to predict ages based on length (or sizes, more generally). The
#' model types that are currently available are 'age-length key', a 'smart
#' age-length key' (see Details), random forest, and gradient boosting machine.
#'
#' \code{fit_age_model} will take the provided length-at-age data and create a
#' model for predicting ages based on length (or whatever is specific as the
#' \code{size_col} argument). The different methods for doing this are:
#' \itemize{
#'   \item{alk}{ -- a basic age-length key}
#'   \item{smart_alk}{ -- a more advanced method that creates an
#'   age-length key from data when enough is available, but borrows data from
#'   elsewhere when it's not. This option requires that you specify levels at
#'   which to fit each alk (i.e. species, county, waterbody, year, etc.).
#'   As an example, if species, county, county, and waterbody are each specified
#'   as levels, then \code{fit_age_model} will create age-length keys for each
#'   waterbody, for each county, and ultimately for each species. The highest
#'   level most general ALK is created at the first level specified.}
#'   \item{rf}{ -- a random forest model that includes any levels specified in
#'   \code{levels}}
#'   \item{gbm}{ -- a gradient boosting machine model that includes any levels
#'   specified in \code{levels}}
#' }
#'
#' For the rf and gbm models, the strings specified in levels get converted to
#' and added to the formula passed those models. For example, if
#' levels = c("spp", "county", "waterbody"), then the resulting formula in the
#' call to the model would be age ~ length + spp + county + waterbody.
#'
#' @param data An object of class data.frame
#' @param model Character. The model type to fit. Options are 'alk',
#' 'halk', 'rf' (for random forest), 'gbm' (for gradient boosting machine)
#' @param levels Character vector. The levels that the age estimating model
#' will fit to. Each level specified must correspond to a column in \code{data}.
#' For models 'alk' and 'halk' an age-length key gets created
#' at each level. For 'rf' and 'gbm' the levels get converted to a formula
#' (see Details)
#' @param age_col Character. The name of the column in data that contains ages
#' @param size_col Character. The name of the column in data that contains sizes
#' @param ... Additional arguments passed onto the various methods
#'
#' @return An object of the appropriate model type according to what was
#' provided by the \code{model} argument
#' @export
#'
#' @examples
#' spp_halk <- fit_age_model(spp_data, levels = "spp")
fit_age_model <- function(data,
                          model = "halk",
                          levels = NULL,
                          age_col = "age",
                          size_col = "length",
                          ...) {
  UseMethod("fit_age_model", data)
}

#' @export
fit_age_model.data.frame <- function(data,
                                     model = "halk",
                                     levels = NULL,
                                     age_col = "age",
                                     size_col = "length",
                                     ...) {
  stopifnot(is.data.frame(data))
  check_length_data(data, size_col)
  check_age_data(data, age_col)
  levels <- add_spp_level(data, levels)
  model <- check_model_type(model)
  if (model == "halk") {
    if (is.null(levels)) {
      stop(
        "If you are trying to create an HALK you need to specify something",
        " in the levels argument."
      )
    } else {
      model <- "alk"
    }
  }
  class(data) <- c(model, "data.frame")
  return(fit_age_model(
    data, levels = levels,
    size_col = size_col,
    age_col = age_col,
    ...)
  )
}

#' @export
fit_age_model.alk <- function(data,
                              model,
                              levels = NULL,
                              age_col = age_col,
                              size_col = size_col,
                              ...) {
  if (is.null(levels)) {
    out <- make_alk(data, sizecol = size_col, agecol = age_col, ...)
    if (is.null(out)) {
      return(NULL)
    } else {
      out <- tibble::as_tibble(out)
      class(out) <- c("alk", class(out))
      return(out)
    }
  } else {
    if (any(!(levels %in% names(data)))) {
      stop("Each level supplied must match a column in the data.")
    }
    temp_levels <- levels
    alks <-
      tibble::as_tibble(data[FALSE, ]) %>%
      dplyr::select(tidyselect::all_of(levels))
    for (i in rev(seq_along(temp_levels))) {
      grouping <- rlang::syms(temp_levels)
      # level_alks <-
      #   data %>%
      #   adjust_plus_min_ages_df(minage = min_age, pls_grp = plus_group) %>%
      #   min_count_laa_data(
      #     temp_levels,
      #     min_age_sample_size,
      #     min_total_sample_size
      #   ) %>%
      #   min_age_groups(temp_levels, min_age_groups)
      # if (is.null(level_alks)) {
      #   next
      # } else {
        level_alks <-
          data %>%
          dplyr::group_by(!!!grouping) %>%
          tidyr::nest() %>%
          dplyr::mutate(alk = purrr::map(
            .data$data,
            make_alk,
            sizecol = size_col,
            agecol = age_col,
            warnings = FALSE,
            ...
          ))
        alks <- dplyr::bind_rows(alks, level_alks)
        temp_levels <- temp_levels[-i]
      # }
    }
    alks <-
      alks %>%
      dplyr::filter(!purrr::map_lgl(.data$alk, is.null))
    if (nrow(alks) == 0) {
      warning(
        "Your HALK did not compute. You probably did not have enough \n  ",
        "samples in a particular age group. Either adjust \n  ",
        "min_sample_size or set min_age or plus_group."
      )
      return(NULL)
    } else {
      out <-
        dplyr::select(alks, tidyselect::all_of(levels), "alk") %>%
        dplyr::arrange(!!!rlang::syms(levels)) %>%
        tibble::as_tibble()
      attr(out, "levels") <- levels
      attr(out, "size_col") <- size_col
      attr(out, "age_col") <- age_col
      # attr(out, "autobin") <- autobin
      # attr(out, "size_bin") <- size_bin
      class(out) <- c("halk_fit", class(out))
      return(out)
    }
  }
}

#' Create a hierarchical age-length key (HALK)
#'
#' This function creates a hierarchically nested age-length key that can be
#' used to estimate age of an organism based on proportion of sampled organisms
#' in each age group.
#'
#' @param data A data.frame with age and size samples
#' @param levels Character vector specifying the levels for HALK creation
#' @param age_col Optional. String of the column name in \code{data} housing
#' age data
#' @param size_col Optional. String of the column name in \code{data} housing
#' size data
#' @param ... Additional arguments passed to \code{\link{make_alk}}
#'
#' @return A \code{\link[tibble]{tibble}} with columns for each level and
#' a column called alk that houses the age-length key for that particular level
#' @export
#'
#' @examples
#' make_halk(spp_data, levels = "spp")
make_halk <- function(data,
                      levels,
                      age_col = "age",
                      size_col = "length",
                      ...) {
  if (is.null(levels)) {
    stop("You must provide levels to create a HALK.")
  } else {
    if (any(!(levels %in% names(data)))) {
      stop("Each level supplied must match a column in the data.")
    }
    temp_levels <- levels
    alks <-
      tibble::as_tibble(data[FALSE, ]) %>%
      dplyr::select(tidyselect::all_of(levels))
    for (i in rev(seq_along(temp_levels))) {
      grouping <- rlang::syms(temp_levels)
      level_alks <-
        data %>%
        dplyr::group_by(!!!grouping) %>%
        tidyr::nest() %>%
        dplyr::mutate(alk = purrr::map(
          .data$data,
          make_alk,
          sizecol = size_col,
          agecol = age_col,
          warnings = FALSE,
          ...
        ))
      alks <- dplyr::bind_rows(alks, level_alks)
      temp_levels <- temp_levels[-i]
    }
    alks <-
      alks %>%
      dplyr::filter(!purrr::map_lgl(.data$alk, is.null))
    if (nrow(alks) == 0) {
      warning(
        "Your HALK did not compute. You probably did not have enough \n  ",
        "samples in a particular age group. Either adjust \n  ",
        "min_sample_size or set min_age or plus_group."
      )
      return(NULL)
    } else {
      out <-
        dplyr::select(alks, tidyselect::all_of(levels), "alk") %>%
        dplyr::arrange(!!!rlang::syms(levels)) %>%
        tibble::as_tibble()
      attr(out, "levels") <- levels
      attr(out, "size_col") <- size_col
      attr(out, "age_col") <- age_col
      class(out) <- c("halk_fit", class(out))
      return(out)
    }
  }
}

