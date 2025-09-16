# modules/mod_armpriors.R
# Prior & posterior Beta densities per arm (faceted)

mod_armpriors_ui <- function(id, height = 320) {
  ns <- NS(id)
  tagList(
    h4("Priors and posteriors by arm"),
    plotOutput(ns("priors_plot"), height = height)
  )
}

mod_armpriors_server <- function(id,
                                 p_t, n_t,       # reactives: treatment truth in [0,1], n_t
                                 p_c, n_c,       # reactives: control truth in [0,1], n_c
                                 prior,          # reactive: "flat" or "power"
                                 prior_args) {   # reactive: list(a_base,b_base,a0,y_0,n_0)

  moduleServer(id, function(input, output, session) {

    # small helper
    `%||%` <- function(x, y) if (is.null(x)) y else x

    # compute Beta params for priors & "expected" posteriors
    params <- reactive({
      pa <- prior_args() %||% list()
      pr <- prior()

      # Baseline priors
      a_base <- pa$a_base %||% 1
      b_base <- pa$b_base %||% 1

      # Treatment arm prior: baseline only (adjust if you later expose a_t/a_b for treatment)
      a_t_prior <- 1
      b_t_prior <- 1

      # Control arm prior: baseline, plus borrowed history if power prior
      a0 <- (pa$a0 %||% 0)
      y0 <- (pa$y_0 %||% 0)
      n0 <- (pa$n_0 %||% 0)

      if (identical(pr, "power")) {
        a_c_prior <- a_base + a0 * y0
        b_c_prior <- b_base + a0 * (n0 - y0)
      } else {
        a_c_prior <- a_base
        b_c_prior <- b_base
      }

      # Expected posteriors based on current inputs (ŷ = round(n * p))
      nt <- as.integer(n_t())
      nc <- as.integer(n_c())
      pt <- p_t()
      pc <- p_c()

      y_t_hat <- round(nt * pt)
      y_c_hat <- round(nc * pc)

      a_t_post <- a_t_prior + y_t_hat
      b_t_post <- b_t_prior + (nt - y_t_hat)

      a_c_post <- a_c_prior + y_c_hat
      b_c_post <- b_c_prior + (nc - y_c_hat)

      list(
        tr_prior = c(a = a_t_prior, b = b_t_prior),
        tr_post  = c(a = a_t_post,  b = b_t_post),
        ct_prior = c(a = a_c_prior, b = b_c_prior),
        ct_post  = c(a = a_c_post,  b = b_c_post),
        meta = list(
          pt = pt, pc = pc, nt = nt, nc = nc,
          prior = pr, a0 = a0, y0 = y0, n0 = n0,
          a_base = a_base, b_base = b_base
        )
      )
    })

    output$priors_plot <- renderPlot({
      par <- params()

      # density grid
      x <- seq(1e-4, 1 - 1e-4, length.out = 501)

      # build tidy data for ggplot
      df <- rbind(
        data.frame(
          arm  = "Treatment",
          dist = "Prior",
          x = x,
          density = dbeta(x, par$tr_prior["a"], par$tr_prior["b"])
        ),
        data.frame(
          arm  = "Treatment",
          dist = "Posterior",
          x = x,
          density = dbeta(x, par$tr_post["a"], par$tr_post["b"])
        ),
        data.frame(
          arm  = "Control",
          dist = "Prior",
          x = x,
          density = dbeta(x, par$ct_prior["a"], par$ct_prior["b"])
        ),
        data.frame(
          arm  = "Control",
          dist = "Posterior",
          x = x,
          density = dbeta(x, par$ct_post["a"], par$ct_post["b"])
        )
      )

      df$arm  <- factor(df$arm,  levels = c("Control", "Treatment"))
      df$dist <- factor(df$dist, levels = c("Prior", "Posterior"))

      # vertical reference lines at assumed truths
      vlines <- data.frame(
        arm = c("Control", "Treatment"),
        x   = c(par$meta$pc, par$meta$pt),
        lab = c(sprintf("p_c = %.0f%%", 100*par$meta$pc),
                sprintf("p_t = %.0f%%", 100*par$meta$pt))
      )

      ggplot(df, aes(x = x, y = density, color = dist, fill = dist)) +
        geom_line(data = subset(df, dist == "Prior"), linewidth = 0.9) +
        geom_area(data = subset(df, dist == "Posterior"),
                  alpha = 0.3, position = "identity") +
        geom_line(data = subset(df, dist == "Posterior"), linewidth = 1.1) +
        geom_vline(data = vlines, aes(xintercept = x), linetype = "dashed") +
        facet_wrap(~ arm, nrow = 1, scales = "fixed") +
        scale_color_manual(values = c("Prior" = "grey40", "Posterior" = "firebrick")) +
        scale_fill_manual(values  = c("Prior" = NA,       "Posterior" = "firebrick")) +
        labs(
          x = expression(theta),
          y = "Density",
          subtitle = {
            m <- par$meta
            if (identical(m$prior, "power")) {
              sprintf("Control prior: Beta(%.1f, %.1f) + a0·History (a0=%.2f, y0=%d, n0=%d); Treatment prior: Beta(1,1).",
                      m$a_base, m$b_base, m$a0, m$y0, m$n0)
            } else {
              sprintf("Control prior: Beta(%.1f, %.1f); Treatment prior: Beta(1,1).",
                      m$a_base, m$b_base)
            }
          }
        ) +
        guides(color = guide_legend(title = NULL),
               fill  = guide_legend(title = NULL)) +
        theme_minimal(base_size = 12) +
        theme(
          strip.text = element_text(face = "bold"),
          legend.position = "bottom"
        )
    })
  })
}
