#' Cortisol ELISA 4PL Standard Curve
#'
#' 输入 6 个标准品原始 OD450，先转换为 OD%，再基于 OD% 拟合 4 参数法（4PL）标准曲线。
#' 适用于吸光度越高、浓度越低的竞争法 ELISA。
#'
#' 固定浓度顺序为：
#' 0, 0.41, 1.10, 4.69, 19.3, 82.8 (nmol/L)
#'
#' @param od450 长度为 6 的数值向量，对应固定浓度顺序下的原始 OD450。
#' @param auto_print 是否自动打印结果，默认 TRUE。
#' @param auto_plot 是否自动显示图形，默认 TRUE。
#'
#' @return 返回一个列表，包括数据、参数、R2、100%OD、Excel公式、图和模型。
#' @export
cortisol1 <- function(od450, auto_print = TRUE, auto_plot = TRUE) {
  if (!is.numeric(od450)) {
    stop("od450 必须是数值向量。")
  }

  if (length(od450) != 6) {
    stop("od450 必须正好包含 6 个数值，顺序对应浓度 0, 0.41, 1.10, 4.69, 19.3, 82.8。")
  }

  if (any(is.na(od450))) {
    stop("od450 不能包含 NA。")
  }

  if (od450[1] <= 0) {
    stop("第1个OD值必须大于0，因为它作为100%OD。")
  }

  conc <- c(0, 0.41, 1.10, 4.69, 19.3, 82.8)
  od100 <- od450[1]

  # 原始OD -> OD%
  od_percent <- (od450 / od100) * 100

  # 对数坐标轴不能包含0，因此仅在拟合和作图时把0替换为 0.01
  conc_fit <- conc
  conc_fit[conc_fit == 0] <- 0.01

  dat <- data.frame(
    conc = conc,
    conc_fit = conc_fit,
    od450 = od450,
    od_percent = od_percent
  )

  # LL.4 参数含义：b=slope, c=lower, d=upper, e=ED50
  # 映射到公式 y = D + (A - D) / (1 + (x / C)^B)
  # B=b, D=c, A=d, C=e
  fit <- drc::drm(
    od_percent ~ conc_fit,
    data = dat,
    fct = drc::LL.4(names = c("B", "D", "A", "C"))
  )

  cf <- stats::coef(fit)
  B <- unname(cf[1])
  D <- unname(cf[2])
  A <- unname(cf[3])
  C <- unname(cf[4])

  pred <- stats::predict(fit, newdata = data.frame(conc_fit = conc_fit))
  sse <- sum((od_percent - pred)^2)
  sst <- sum((od_percent - mean(od_percent))^2)
  r2 <- 1 - sse / sst

  equation <- "y = D + (A - D) / (1 + (x / C)^B)"
  inverse_formula_percent <- "x = C * (((A - D) / (y - D) - 1)^(1 / B))"
  inverse_formula_raw <- paste0(
    "x = C * (((A - D) / ((((OD/", signif(od100, 10), ")*100) - D)) - 1)^(1 / B))"
  )

  num_dot <- function(x) format(x, digits = 15, scientific = FALSE, trim = TRUE)
  num_comma <- function(x) chartr(".", ",", num_dot(x))

  # Excel 原始OD公式：A1 填原始OD
  excel_formula_raw_de <- sprintf(
    '=WENNFEHLER(%s*(((%s-%s)/(((A1/%s)*100)-%s)-1)^(1/%s));"")',
    num_comma(C), num_comma(A), num_comma(D), num_comma(od100), num_comma(D), num_comma(B)
  )

  excel_formula_raw_cn <- sprintf(
    '=IFERROR(%s*(((%s-%s)/(((A1/%s)*100)-%s)-1)^(1/%s)),"")',
    num_dot(C), num_dot(A), num_dot(D), num_dot(od100), num_dot(D), num_dot(B)
  )

  # Excel 百分比OD公式：A1 填 OD%
  excel_formula_percent_de <- sprintf(
    '=WENNFEHLER(%s*(((%s-%s)/(A1-%s)-1)^(1/%s));"")',
    num_comma(C), num_comma(A), num_comma(D), num_comma(D), num_comma(B)
  )

  excel_formula_percent_cn <- sprintf(
    '=IFERROR(%s*(((%s-%s)/(A1-%s)-1)^(1/%s)),"")',
    num_dot(C), num_dot(A), num_dot(D), num_dot(D), num_dot(B)
  )

  predict_conc <- function(sample_od) {
    if (!is.numeric(sample_od)) {
      stop("sample_od 必须是数值。")
    }

    sample_percent <- (sample_od / od100) * 100
    x <- C * (((A - D) / (sample_percent - D) - 1)^(1 / B))
    x[!is.finite(x)] <- NA_real_
    x[x < 0] <- NA_real_
    x
  }

  xgrid <- exp(seq(log(0.01), log(max(conc_fit)), length.out = 300))
  ygrid <- stats::predict(fit, newdata = data.frame(conc_fit = xgrid))
  plot_df <- data.frame(conc_fit = xgrid, od_percent = ygrid)
  plot_df <- plot_df[order(plot_df$conc_fit), ]

  formula_expr <- bquote(y == D + frac(A - D, 1 + (x / C)^B))
  r2_expr <- bquote(R^2 == .(signif(r2, 5)))

  other_text <- paste0(
    "A = ", signif(A, 6), "\n",
    "B = ", signif(B, 6), "\n",
    "C = ", signif(C, 6), "\n",
    "D = ", signif(D, 6), "\n",
    "100%OD = ", signif(od100, 6), "\n",
    "OD%反算: x = C * (((A - D) / (y - D) - 1)^(1 / B))\n",
    "原始OD反算: x = C * (((A - D) / ((((OD/", signif(od100, 6), ")*100) - D)) - 1)^(1 / B))"
  )

  y_min <- min(c(dat$od_percent, plot_df$od_percent), na.rm = TRUE)
  y_max <- max(c(dat$od_percent, plot_df$od_percent), na.rm = TRUE)
  y_range <- y_max - y_min

  x_left <- 0.012
  x_right <- max(plot_df$conc_fit, na.rm = TRUE) * 0.98

  y_formula <- y_min + 0.95 * y_range
  y_r2 <- y_min + 0.84 * y_range
  y_other <- y_min + 0.02 * y_range

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = conc_fit, y = od_percent)) +
    ggplot2::geom_point(size = 3, color = "darkblue") +
    ggplot2::geom_line(
      data = plot_df,
      ggplot2::aes(x = conc_fit, y = od_percent),
      linewidth = 1,
      color = "darkblue"
    ) +
    ggplot2::scale_x_log10(
      breaks = c(0.01, 0.1, 1, 10, 100),
      labels = c("0.01", "0.1", "1", "10", "100")
    ) +
    ggplot2::annotate(
      "text",
      x = x_right,
      y = y_formula,
      label = deparse(formula_expr),
      parse = TRUE,
      hjust = 1,
      vjust = 1,
      size = 5,
      color = "darkred"
    ) +
    ggplot2::annotate(
      "text",
      x = x_right,
      y = y_r2,
      label = deparse(r2_expr),
      parse = TRUE,
      hjust = 1,
      vjust = 1,
      size = 5,
      color = "darkred"
    ) +
    ggplot2::annotate(
      "text",
      x = x_left,
      y = y_other,
      label = other_text,
      hjust = 0,
      vjust = 0,
      size = 3.2,
      color = "black"
    ) +
    ggplot2::labs(
      title = "Cortisol ELISA Standard Curve (4PL)",
      x = "Concentration (nmol/L, log scale)",
      y = "OD%"
    ) +
    ggplot2::theme_bw()

  res <- list(
    data = dat[, c("conc", "od450", "od_percent")],
    parameters = c(A = A, B = B, C = C, D = D),
    r2 = r2,
    od100 = od100,
    equation = equation,
    inverse_formula_percent = inverse_formula_percent,
    inverse_formula_raw = inverse_formula_raw,
    excel_formula_raw = list(
      de = excel_formula_raw_de,
      cn = excel_formula_raw_cn
    ),
    excel_formula_percent = list(
      de = excel_formula_percent_de,
      cn = excel_formula_percent_cn
    ),
    predict_conc = predict_conc,
    plot = p,
    model = fit
  )

  if (auto_print) {
    cat("\n===== ggELISA::cortisol1 结果 =====\n")
    cat("A = ", signif(A, 6), "\n", sep = "")
    cat("B = ", signif(B, 6), "\n", sep = "")
    cat("C = ", signif(C, 6), "\n", sep = "")
    cat("D = ", signif(D, 6), "\n", sep = "")
    cat("R² = ", signif(r2, 6), "\n", sep = "")
    cat("100%OD = ", signif(od100, 6), "\n\n", sep = "")
    cat("OD%反算公式:\n", inverse_formula_percent, "\n\n", sep = "")
    cat("原始OD反算公式:\n", inverse_formula_raw, "\n\n", sep = "")
    cat("Excel原始OD公式（德版）:\n", excel_formula_raw_de, "\n\n", sep = "")
    cat("Excel原始OD公式（中文版）:\n", excel_formula_raw_cn, "\n\n", sep = "")
    cat("Excel百分比OD公式（德版）:\n", excel_formula_percent_de, "\n\n", sep = "")
    cat("Excel百分比OD公式（中文版）:\n", excel_formula_percent_cn, "\n", sep = "")
  }

  if (auto_plot) {
    print(p)
  }

  invisible(res)
}
