#!/usr/bin/env Rscript

# Deterministic Nature-style renderer for CMM schema-v2 transformation runs (SC-02).
# Scientific decisions are made upstream. This script only reads tidy artifacts declared in
# 00_manifest.json and maps their existing columns to publication graphics.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: render_transformation_figures.R RUN_DIR MANIFEST OUTPUT_DIR FIGURE_MANIFEST")
}

minimum_r_version <- "4.3.2"
minimum_package_versions <- c(
  jsonlite = "1.8.8",
  ggplot2 = "3.5.2",
  ggrepel = "0.9.6",
  svglite = "2.2.1",
  ragg = "1.2.7"
)
required_packages <- names(minimum_package_versions)
if (getRversion() < numeric_version(minimum_r_version)) {
  stop(sprintf(
    paste0(
      "incompatible R version: actual %s; required >= %s. ",
      "Restore the checked-in renv.lock or install a compatible R runtime."
    ),
    as.character(getRversion()),
    minimum_r_version
  ))
}
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop(sprintf("missing required R packages: %s", paste(missing_packages, collapse = ", ")))
}
incompatible_packages <- required_packages[vapply(
  required_packages,
  function(package) utils::packageVersion(package) < numeric_version(minimum_package_versions[[package]]),
  logical(1)
)]
if (length(incompatible_packages) > 0L) {
  details <- vapply(
    incompatible_packages,
    function(package) sprintf(
      "%s actual %s; required >= %s",
      package,
      as.character(utils::packageVersion(package)),
      minimum_package_versions[[package]]
    ),
    character(1)
  )
  stop(sprintf(
    "incompatible R package versions: %s. Restore the checked-in renv.lock.",
    paste(details, collapse = "; ")
  ))
}
# A publication render must not silently discard rows or use an invalid scale/device.
options(warn = 2)

run_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
manifest_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[3]], winslash = "/", mustWork = FALSE)
figure_manifest_path <- normalizePath(args[[4]], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(left, right) {
  if (is.null(left) || length(left) == 0L) right else left
}

manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
if (!identical(manifest$schema_version, 2L)) {
  stop("manifest schema_version must be 2")
}
if (!identical(as.character(manifest$workflow %||% ""), "transformation_target_discovery")) {
  stop("manifest is not a transformation_target_discovery run")
}

artifact_entry <- function(role) {
  manifest$artifacts[[role]]
}

artifact_status <- function(role) {
  entry <- artifact_entry(role)
  if (is.null(entry)) return("missing")
  as.character(entry$status %||% "missing")
}

artifact_relative <- function(role) {
  entry <- artifact_entry(role)
  if (is.null(entry)) return(NULL)
  path <- entry$path %||% NULL
  if (is.null(path) || !nzchar(as.character(path))) return(NULL)
  as.character(path)
}

artifact_available <- function(role) {
  artifact_status(role) %in% c("complete", "partial") && !is.null(artifact_relative(role))
}

artifact_path <- function(role) {
  relative <- artifact_relative(role)
  if (is.null(relative)) stop(sprintf("artifact role %s has no path", role))
  file.path(run_dir, relative)
}

read_artifact_csv <- function(role) {
  utils::read.csv(
    artifact_path(role),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN")
  )
}

relative_to_run <- function(path) {
  absolute <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(run_dir, "/")
  if (!startsWith(absolute, prefix)) stop(sprintf("output escaped run directory: %s", absolute))
  substring(absolute, nchar(prefix) + 1L)
}

# Which candidate to mark is a reading choice, not a property of the run, so it arrives the same
# way the font does rather than being written into the bundle.
highlight <- Sys.getenv("CMM_TRANSFORMATION_HIGHLIGHT", unset = "")
highlight <- if (nzchar(trimws(highlight))) trimws(highlight) else NULL

font_requested <- Sys.getenv("CMM_FIGURE_FONT", unset = "Helvetica")
font_family <- font_requested
# The base PDF device provides editable vector text without a Cairo/X11 runtime dependency.
if (!(font_family %in% names(grDevices::pdfFonts()))) {
  font_family <- "Helvetica"
}
font_resolved <- "not inspected"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  match <- tryCatch(systemfonts::match_fonts(font_requested), error = function(e) NULL)
  if (!is.null(match) && nrow(match) > 0L && nzchar(match$path[[1]])) {
    font_resolved <- basename(match$path[[1]])
  } else {
    font_family <- "Helvetica"
    fallback <- tryCatch(systemfonts::match_fonts(font_family), error = function(e) NULL)
    if (!is.null(fallback) && nrow(fallback) > 0L) {
      font_resolved <- basename(fallback$path[[1]])
    }
  }
}

blue <- "#1769AA"
orange <- "#D55E00"
grey <- "#767676"
light_grey <- "#D9D9D9"

theme_nature <- function() {
  ggplot2::theme_classic(base_family = font_family, base_size = 6) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = grid::unit(1.2, "mm"),
      axis.title = ggplot2::element_text(size = 6.5, colour = "black"),
      axis.text = ggplot2::element_text(size = 6, colour = "black"),
      legend.title = ggplot2::element_text(size = 6),
      legend.text = ggplot2::element_text(size = 5.5),
      legend.key.height = grid::unit(3, "mm"),
      legend.key.width = grid::unit(4, "mm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 6.5, face = "bold", colour = "black"),
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(3, 3, 3, 3, unit = "mm")
    )
}

draw_on_device <- function(open_device, plot) {
  open_device()
  device <- grDevices::dev.cur()
  tryCatch(
    print(plot),
    finally = {
      if (grDevices::dev.cur() == device) grDevices::dev.off()
    }
  )
}

save_triplet <- function(plot, figure_id, width_mm, height_mm) {
  paths <- list(
    png = file.path(output_dir, paste0(figure_id, ".png")),
    pdf = file.path(output_dir, paste0(figure_id, ".pdf")),
    svg = file.path(output_dir, paste0(figure_id, ".svg"))
  )
  draw_on_device(
    function() ragg::agg_png(
      filename = paths$png,
      width = width_mm,
      height = height_mm,
      units = "mm",
      res = 300,
      background = "white"
    ),
    plot
  )
  draw_on_device(
    function() grDevices::pdf(
      file = paths$pdf,
      width = width_mm / 25.4,
      height = height_mm / 25.4,
      family = font_family,
      onefile = TRUE,
      useDingbats = FALSE,
      paper = "special"
    ),
    plot
  )
  draw_on_device(
    function() svglite::svglite(
      file = paths$svg,
      width = width_mm / 25.4,
      height = height_mm / 25.4,
      bg = "white"
    ),
    plot
  )
  lapply(paths, relative_to_run)
}

figures <- list()
required_failed <- FALSE

record_skipped <- function(id, label, caption, roles, reason, required = FALSE) {
  figures[[length(figures) + 1L]] <<- list(
    id = id,
    label = label,
    status = if (required) "failed" else "skipped",
    caption = caption,
    alt = caption,
    sources = Filter(Negate(is.null), lapply(roles, artifact_relative)),
    reason = reason
  )
  if (required) required_failed <<- TRUE
}

record_rendered <- function(id, label, caption, alt, roles, width_mm, height_mm, plot) {
  outputs <- save_triplet(plot, id, width_mm, height_mm)
  figures[[length(figures) + 1L]] <<- list(
    id = id,
    label = label,
    status = "rendered",
    caption = caption,
    alt = alt,
    width_mm = width_mm,
    height_mm = height_mm,
    dpi = 300,
    sources = unname(lapply(roles, artifact_relative)),
    outputs = outputs
  )
}

render_panel <- function(id, label, caption, alt, roles, width_mm, height_mm, required, build) {
  unavailable <- roles[!vapply(roles, artifact_available, logical(1))]
  if (length(unavailable) > 0L) {
    reason <- sprintf("Unavailable artifact roles: %s", paste(unavailable, collapse = ", "))
    record_skipped(id, label, caption, roles, reason, required = required)
    return(invisible(NULL))
  }
  tryCatch(
    {
      plot <- build()
      rendered_caption <- attr(plot, "cmm_caption", exact = TRUE) %||% caption
      rendered_alt <- attr(plot, "cmm_alt", exact = TRUE) %||% alt
      record_rendered(id, label, rendered_caption, rendered_alt, roles, width_mm, height_mm, plot)
    },
    error = function(error) {
      record_skipped(id, label, caption, roles, conditionMessage(error), required = required)
    }
  )
}

highlight_flag <- function(ids) {
  if (is.null(highlight)) rep(FALSE, length(ids)) else ids == highlight
}

# A gene symbol or reaction name where the analysis stage resolved one, and the identifier
# otherwise. Runs written before target_name existed carry no such column and fall back
# cleanly, which is why this reads the frame rather than requiring the column.
display_label <- function(frame) {
  ids <- as.character(frame$target_id)
  if (!"target_name" %in% names(frame)) return(ids)
  names_col <- as.character(frame$target_name)
  names_col[is.na(names_col)] <- ""
  ifelse(nzchar(trimws(names_col)), trimws(names_col), ids)
}

# --- Figure 1 — best-case against worst-case component score ----------------------------
render_panel(
  id = "fig01_component_scores",
  label = "Figure 1",
  caption = paste(
    "Best-case against worst-case transformation score for every scored candidate.",
    "Candidates in the lower-right quadrant move flux toward the target state even in the",
    "worst case, and are the only ones the robust score amplifies."
  ),
  alt = "Scatter of best-case against worst-case transformation score, divided into quadrants",
  roles = c("transformation_ranking"),
  width_mm = 89,
  height_mm = 65,
  # Not required: under MTA the three component scores are one number, so this panel has
  # nothing to draw, and a run that legitimately cannot produce it is still a complete run.
  required = FALSE,
  build = function() {
    frame <- read_artifact_csv("transformation_ranking")
    needed <- c("bTS", "mTS", "wTS")
    if (!all(needed %in% names(frame))) {
      stop("the ranking carries no bTS/mTS/wTS columns, which only rMTA reports")
    }
    for (col in needed) frame[[col]] <- as.numeric(frame[[col]])
    frame <- frame[stats::complete.cases(frame[, needed, drop = FALSE]), , drop = FALSE]
    if (nrow(frame) == 0L) stop("the ranking contains no scored candidate")
    if (isTRUE(all.equal(frame$bTS, frame$wTS))) {
      stop(paste(
        "best-case and worst-case scores are identical, as they are under MTA;",
        "the contrast is defined only for rMTA"
      ))
    }
    frame$amplified <- frame$bTS > 0 & frame$mTS > 0 & frame$wTS < 0
    frame$label <- display_label(frame)
    marked <- frame[frame$amplified | highlight_flag(as.character(frame$target_id)), , drop = FALSE]
    legend_labels <- c(`FALSE` = "Robust score = MOMA component", `TRUE` = "Amplified by the robust score")

    plot <- ggplot2::ggplot(frame, ggplot2::aes(x = .data$bTS, y = .data$wTS)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = grey, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = grey, linetype = "dashed") +
      # Amplified candidates are few and sit close to wTS = 0, where the bulk of the cloud
      # also crowds. Drawing them as a second, larger layer keeps them visible whatever the
      # axis range turns out to be.
      ggplot2::geom_point(
        data = frame[!frame$amplified, , drop = FALSE],
        ggplot2::aes(colour = .data$amplified, shape = .data$amplified), size = 0.8
      ) +
      ggplot2::geom_point(
        data = frame[frame$amplified, , drop = FALSE],
        ggplot2::aes(colour = .data$amplified, shape = .data$amplified), size = 1.8
      ) +
      ggplot2::scale_colour_manual(
        values = c(`FALSE` = blue, `TRUE` = orange), labels = legend_labels,
        breaks = c(FALSE, TRUE), drop = FALSE, name = NULL
      ) +
      ggplot2::scale_shape_manual(
        values = c(`FALSE` = 16, `TRUE` = 17), labels = legend_labels,
        breaks = c(FALSE, TRUE), drop = FALSE, name = NULL
      ) +
      ggplot2::labs(x = "Best-case score (bTS)", y = "Worst-case score (wTS)") +
      theme_nature() +
      ggplot2::theme(legend.position = "bottom")
    if (nrow(marked) > 0L) {
      plot <- plot +
        ggrepel::geom_text_repel(
          data = marked,
          mapping = ggplot2::aes(label = .data$label),
          size = 2,
          colour = "black",
          family = font_family,
          min.segment.length = 0,
          segment.size = 0.25,
          max.overlaps = 20L,
          seed = 0L
        )
    }
    caption <- sprintf(
      paste(
        "Best-case (bTS) against worst-case (wTS) transformation score for %d scored candidates.",
        "The %d in the lower-right quadrant move flux toward the target state even in the worst",
        "case and are the only ones the robust score amplifies; every other candidate carries",
        "its MOMA-based component unchanged."
      ),
      nrow(frame),
      sum(frame$amplified)
    )
    attr(plot, "cmm_caption") <- caption
    attr(plot, "cmm_alt") <- caption
    plot
  }
)

# --- Figure 2 — transformation score against rank ---------------------------------------
render_panel(
  id = "fig02_transformation_ranking",
  label = "Figure 2",
  caption = paste(
    "Transformation score against rank. A score that decays smoothly rather than falling away",
    "means the cut between the leading candidates and the rest is a choice, not a boundary."
  ),
  alt = "Scatter of transformation score against rank for every scored candidate",
  roles = c("transformation_ranking"),
  width_mm = 89,
  height_mm = 60,
  required = TRUE,
  build = function() {
    frame <- read_artifact_csv("transformation_ranking")
    frame <- frame[!is.na(frame$score) & !is.na(frame$rank), , drop = FALSE]
    if (nrow(frame) == 0L) stop("the ranking contains no scored candidate")
    frame$rank <- as.integer(frame$rank)
    frame$score <- as.numeric(frame$score)
    frame$label <- display_label(frame)
    frame$marked <- highlight_flag(as.character(frame$target_id))
    marked <- frame[frame$marked, , drop = FALSE]

    plot <- ggplot2::ggplot(frame, ggplot2::aes(x = .data$rank, y = .data$score)) +
      ggplot2::geom_line(linewidth = 0.3, colour = light_grey) +
      ggplot2::geom_point(size = 0.7, colour = blue) +
      ggplot2::labs(x = "Rank", y = "Transformation score") +
      theme_nature()
    if (nrow(marked) > 0L) {
      plot <- plot +
        ggplot2::geom_point(data = marked, size = 1.4, colour = orange) +
        ggrepel::geom_text_repel(
          data = marked,
          mapping = ggplot2::aes(label = .data$label),
          size = 2,
          colour = orange,
          family = font_family,
          min.segment.length = 0,
          segment.size = 0.25,
          seed = 0L
        )
    }
    caption <- sprintf(
      paste(
        "Transformation score against rank for %d scored %s.",
        "A score that decays smoothly rather than falling away means the cut between the",
        "leading candidates and the rest is a choice, not a boundary."
      ),
      nrow(frame),
      if (nrow(frame) == 1L) "candidate" else "candidates"
    )
    attr(plot, "cmm_caption") <- caption
    attr(plot, "cmm_alt") <- caption
    plot
  }
)

# --- Figure 3 — transformation rank against the MOMA baseline ---------------------------
render_panel(
  id = "fig03_ranking_vs_moma",
  label = "Figure 3",
  caption = paste(
    "Each candidate's rank under the transformation method and under the MOMA baseline.",
    "Points on the diagonal are candidates the two methods agree on."
  ),
  alt = "Scatter of transformation rank against MOMA rank, with the identity line drawn",
  roles = c("transformation_ranking", "moma_baseline"),
  width_mm = 89,
  height_mm = 60,
  required = FALSE,
  build = function() {
    ranking <- read_artifact_csv("transformation_ranking")
    baseline <- read_artifact_csv("moma_baseline")
    scored <- is.finite(as.numeric(baseline$moma_score))
    if ("status" %in% names(baseline)) {
      scored <- scored & !is.na(baseline$status) & baseline$status == "optimal"
    }
    excluded <- sum(!scored)
    baseline <- baseline[scored, , drop = FALSE]
    merged <- merge(
      data.frame(
        target_id = as.character(ranking$target_id),
        label = display_label(ranking),
        transformation_rank = as.integer(ranking$rank),
        stringsAsFactors = FALSE
      ),
      data.frame(
        target_id = as.character(baseline$target_id),
        moma_rank = as.integer(baseline$rank),
        stringsAsFactors = FALSE
      ),
      by = "target_id"
    )
    merged <- merged[stats::complete.cases(merged), , drop = FALSE]
    if (nrow(merged) == 0L) stop("no candidate has a successful MOMA solve for rank comparison")
    merged$marked <- highlight_flag(merged$target_id)
    marked <- merged[merged$marked, , drop = FALSE]
    limit <- max(merged$transformation_rank, merged$moma_rank)

    plot <- ggplot2::ggplot(
      merged,
      ggplot2::aes(x = .data$transformation_rank, y = .data$moma_rank)
    ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = grey, linetype = "dashed") +
      ggplot2::geom_point(size = 0.9, colour = blue) +
      ggplot2::coord_equal(xlim = c(1, limit), ylim = c(1, limit)) +
      ggplot2::labs(x = "Transformation rank", y = "MOMA rank") +
      theme_nature()
    if (nrow(marked) > 0L) {
      plot <- plot +
        ggplot2::geom_point(data = marked, size = 1.5, colour = orange) +
        ggrepel::geom_text_repel(
          data = marked,
          mapping = ggplot2::aes(label = .data$label),
          size = 2,
          colour = orange,
          family = font_family,
          min.segment.length = 0,
          segment.size = 0.25,
          seed = 0L
        )
    }
    caption <- sprintf(
      paste(
        "Ranks for %d candidates with successful MOMA solves; points on the diagonal",
        "are candidates the two methods agree on.",
        "%d unsuccessful MOMA solves are excluded from this comparison and retained in the CSV."
      ),
      nrow(merged), excluded
    )
    attr(plot, "cmm_caption") <- caption
    attr(plot, "cmm_alt") <- caption
    plot
  }
)

# --- Figure 4 — rank against epsilon ------------------------------------------------------
render_panel(
  id = "fig04_epsilon_sensitivity",
  label = "Figure 4",
  caption = paste(
    "How each leading candidate's rank moves with epsilon. A candidate whose rank is flat",
    "across the sweep was not produced by the value that was chosen."
  ),
  alt = "Rank against epsilon for the leading candidates, one line per candidate",
  roles = c("epsilon_sensitivity"),
  width_mm = 89,
  height_mm = 60,
  required = FALSE,
  build = function() {
    frame <- read_artifact_csv("epsilon_sensitivity")
    frame <- frame[!is.na(frame$epsilon) & !is.na(frame$rank), , drop = FALSE]
    if (nrow(frame) == 0L) stop("the epsilon sweep contains no ranked candidate")
    frame$epsilon <- as.numeric(frame$epsilon)
    frame$rank <- as.integer(frame$rank)
    frame$target_id <- as.character(frame$target_id)

    # One line per candidate is unreadable at genome scale; the sweep exists to show whether the
    # leading candidates hold their place, so it is drawn for those and the cut is stated.
    leading <- unique(frame$target_id[frame$rank <= 5L])
    if (!is.null(highlight) && highlight %in% frame$target_id) {
      leading <- unique(c(leading, highlight))
    }
    drawn <- frame[frame$target_id %in% leading, , drop = FALSE]
    if (nrow(drawn) == 0L) stop("no candidate reaches the leading ranks at any swept epsilon")

    plot <- ggplot2::ggplot(
      drawn,
      ggplot2::aes(
        x = .data$epsilon,
        y = .data$rank,
        colour = .data$target_id,
        group = .data$target_id
      )
    ) +
      ggplot2::geom_line(linewidth = 0.4) +
      ggplot2::geom_point(size = 0.8) +
      ggplot2::scale_y_reverse() +
      ggplot2::scale_colour_manual(values = grDevices::hcl.colors(length(leading), "Dark 3")) +
      ggplot2::labs(x = "Epsilon (flux units)", y = "Rank", colour = NULL) +
      theme_nature()
    caption <- sprintf(
      paste(
        "How rank moves with epsilon for the %d %s reaching the leading five at any swept",
        "value. A candidate whose rank is flat across the sweep was not produced by the value",
        "that was chosen."
      ),
      length(leading),
      if (length(leading) == 1L) "candidate" else "candidates"
    )
    attr(plot, "cmm_caption") <- caption
    attr(plot, "cmm_alt") <- caption
    plot
  }
)

renderer_versions <- lapply(required_packages, function(package) {
  as.character(utils::packageVersion(package))
})
names(renderer_versions) <- required_packages

figure_manifest <- list(
  schema_version = 2L,
  renderer = list(
    engine = "R/ggplot2",
    r = R.version.string,
    minimum_r = minimum_r_version,
    packages = renderer_versions,
    minimum_packages = as.list(minimum_package_versions),
    script_sha256 = Sys.getenv("CMM_RENDERER_SHA256", unset = "unknown"),
    font_requested = font_requested,
    font_family = font_family,
    font_resolved = font_resolved,
    specification = list(
      widths_mm = c(89L, 180L),
      raster_dpi = 300L,
      vector_formats = c("pdf", "svg"),
      colour_space = "RGB"
    )
  ),
  figures = figures
)

jsonlite::write_json(
  figure_manifest,
  path = figure_manifest_path,
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null"
)

if (required_failed) {
  stop("one or more required transformation figures failed; see figures/figure_manifest.json")
}
