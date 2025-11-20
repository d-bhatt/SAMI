library(shiny)
library(ggplot2)
library(patchwork)
library(viridis)
library(dplyr)
library(grid)

theme_set(theme_minimal(base_size = 12, base_family = "Arial"))

# ===========================================================
# --- Helper Functions ---
# ===========================================================

# Initialize spatial positions of cells
initialize_positions <- function(n_cells, geometry = "monolayer", L = 100) {
  if (geometry == "monolayer") {
    x <- runif(n_cells, 0, L)
    y <- runif(n_cells, 0, L)
    return(data.frame(cell_id = 1:n_cells, x = x, y = y))
  } else if (geometry == "spheroid") {
    coords <- matrix(NA, nrow = n_cells, ncol = 3)
    i <- 1
    while (i <= n_cells) {
      point <- runif(3, -L, L)
      if (sum(point^2) <= L^2) {
        coords[i,] <- point
        i <- i + 1
      }
    }
    return(data.frame(cell_id = 1:n_cells, x = coords[,1], y = coords[,2], z = coords[,3]))
  }
}

# Compute crowding factor (neighbor density effects)
compute_crowding <- function(positions, radius = 10,
                             kernel = c("flat","gaussian"),
                             mode   = c("inhibit","promote"),
                             beta   = 1.0,
                             clamp  = c(0.2, 2.0)) {
  kernel <- match.arg(kernel); mode <- match.arg(mode)
  n <- nrow(positions)
  dens <- numeric(n)
  
  for (i in 1:n) {
    dx <- positions$x - positions$x[i]
    dy <- positions$y - positions$y[i]
    if ("z" %in% names(positions)) {
      dz <- positions$z - positions$z[i]
      d  <- sqrt(dx^2 + dy^2 + dz^2)
    } else {
      d  <- sqrt(dx^2 + dy^2)
    }
    if (kernel == "flat") {
      dens[i] <- sum(d > 0 & d < radius)
    } else {
      sigma <- radius/2
      w <- exp(-0.5 * (d/sigma)^2)
      w[d == 0] <- 0
      dens[i] <- sum(w[d < 3*radius])
    }
  }
  
  dens_ref <- max(median(dens), 1e-6)
  if (mode == "inhibit") {
    f <- 1 / (1 + beta * (dens/dens_ref))
  } else {
    f <- 1 + beta * (dens/dens_ref)
  }
  
  f <- pmax(clamp[1], pmin(clamp[2], f))
  list(factor = f, density = dens, density_ref = dens_ref)
}

# Initialize T cells
initialize_tcells <- function(n_tcells, geometry = "spheroid", L = 100, placement = "random") {
  tpos <- initialize_positions(n_tcells, geometry, L)
  tpos$kills <- 0L
  tpos
}

# Core simulation function
times <- seq(0,24,by=1)
simulate_cell <- function(therapy, dose, cell_id, 
                          x = NA, y = NA, z = NULL,
                          crowding_factor = 1,
                          ischemia_rate = 0.05,
                          ischemia_pattern = "none",
                          resistance_prob = 0.1,
                          spatial = FALSE,
                          infection_window = 1,
                          superinfection_resistance = TRUE,
                          delivery_mode = "random",
                          tcell_positions = NULL,
                          kill_radius = 10,
                          antigen_threshold = 1000,
                          recognition_variability = 0.2,
                          max_kills_per_tcell = 10,
                          enable_tcell_killing = TRUE,
                          enable_ischemia_mod = TRUE,
                          enable_crowding_mod = TRUE,
                          geometry = "monolayer") {
  
  mRNA <- 0; replicase <- 0; Protein <- 0; sub_mRNA <- 0; superinfected <- FALSE
  alive <- TRUE; killers <- 0; kill_time <- NA; killer_id <- NA
  
  if (spatial) {
    distance_from_center <- sqrt(x^2 + y^2 + ifelse(!is.null(z), z^2, 0))
    max_dist <- if (!is.null(z)) sqrt(3)*50 else sqrt(2)*50
    
    if (ischemia_pattern == "center") {
      ischemia_factor <- exp(-ischemia_rate * (max_dist - distance_from_center))
    } else if (ischemia_pattern == "periphery") {
      ischemia_factor <- exp(-ischemia_rate * distance_from_center)
    } else {
      ischemia_factor <- 1
    }
    
    if (delivery_mode == "center") {
      dose_eff <- dose * exp(-distance_from_center / 20)
    } else if (delivery_mode == "random") {
      dose_eff <- dose * runif(1, 0.5, 1.5)
    } else if (delivery_mode == "border") {
      dose_eff <- dose * (distance_from_center / max_dist)
    } else {
      dose_eff <- dose
    }
    
    dose_eff <- dose_eff * crowding_factor
  } else {
    dose_eff <- dose
    ischemia_factor <- 1
  }
  
  is_resistant <- runif(1) < resistance_prob
  resistance_factor <- if (is_resistant) 0.3 else 1.0
  
  Protein_max <- if (is_resistant) rnorm(1, mean = 1500, sd = 300) else rnorm(1, mean = 5000, sd = 1000)
  Protein_max <- max(Protein_max, 500)
  transcription_rate <- if (is_resistant) rnorm(1, mean = therapy$rates["k1"] * 0.2, sd = 0.05) else rnorm(1, mean = therapy$rates["k1"], sd = 0.1)
  transcription_rate <- transcription_rate * crowding_factor * ischemia_factor
  transcription_rate <- max(transcription_rate, 0)
  
  out <- data.frame(time = times, Antigen = NA, Antigen_vaccine = NA, Antigen_innate = NA)
  
  if (!is.null(tcell_positions)) local_tcell_kills <- integer(nrow(tcell_positions))
  
  for (t in times) {
    if (!alive) break
    
    if (t <= infection_window) {
      lambda <- dose_eff
      new_mRNA <- rpois(1, lambda)
      mRNA <- mRNA + new_mRNA
      if (superinfection_resistance && new_mRNA > 0) superinfected <- TRUE
    }
    
    replicase <- replicase + rpois(1, transcription_rate * mRNA)
    sub_mRNA <- sub_mRNA + rpois(1, therapy$rates["k2"] * replicase)
    eff_rate <- sub_mRNA * (1 - Protein / Protein_max)
    eff_rate <- eff_rate * resistance_factor * ischemia_factor * crowding_factor
    Protein <- Protein + rpois(1, pmax(eff_rate,0))
    
    out$Antigen[times == t] <- Protein
    out$Antigen_vaccine[times == t] <- Protein
    out$Antigen_innate[times == t] <- rpois(1, 50)
    
    if (enable_tcell_killing && !is.null(tcell_positions) && spatial) {
      dx <- tcell_positions$x - x
      dy <- tcell_positions$y - y
      dz <- if(!is.null(z) && "z" %in% names(tcell_positions)) tcell_positions$z - z else 0
      dists <- sqrt(dx^2 + dy^2 + dz^2)
      nearby_idx <- which(dists <= kill_radius)
      if (length(nearby_idx) > 0) {
        for (ti in nearby_idx) {
          if (local_tcell_kills[ti] >= max_kills_per_tcell) next
          antigen_signal <- out$Antigen_vaccine[times==t]
          threshold <- antigen_threshold * runif(1, 1-recognition_variability, 1+recognition_variability)
          prob <- pmin(1, antigen_signal/(threshold+1))
          if (enable_ischemia_mod) prob <- prob * ischemia_factor
          if (enable_crowding_mod) prob <- prob * crowding_factor
          if (runif(1) < prob) {
            alive <- FALSE
            killers <- killers + 1
            kill_time <- t
            killer_id <- ti
            local_tcell_kills[ti] <- local_tcell_kills[ti] + 1L
            break
          }
        }
      }
    }
  }
  
  out$cell_id <- cell_id
  out$superinfected <- superinfected
  out$Protein_max <- Protein_max
  out$transcription_rate <- transcription_rate
  out$is_resistant <- is_resistant
  out$dose_eff <- dose_eff
  out$alive <- alive
  out$kill_time <- kill_time
  out$killer_id <- killer_id
  if (spatial) {
    out$x <- x; out$y <- y
    if (!is.null(z)) out$z <- z
    out$crowding_factor <- crowding_factor
    out$ischemia_factor <- ischemia_factor
    out$distance <- distance_from_center
  }
  return(out)
}

# Therapy placeholder
therapies <- list(Replicon=list(rates=c(k1=1.0, k2=4.0), name="Replicon"))

# ===========================================================
# --- UI ---
# ===========================================================
ui <- fluidPage(
  titlePanel("sa-mRNA Immunotherapy model (SAMI)"),
  sidebarLayout(
    sidebarPanel(
      numericInput("n_cells", "Tumor Cells:", 1000, min=50, max=2000, step=50),
      numericInput("n_tcells", "T Cells:", 50, min=0, max=200, step=10),
      sliderInput("dose", "Dose:", min=0, max=1000, value=10),
      selectInput("geometry", "Geometry:", choices=c("monolayer","spheroid"), selected="spheroid"),
      checkboxInput("vary_radius", "Vary Spheroid Radius (Fixed # Cells)", value = TRUE), 
      conditionalPanel ( condition = "input.vary_radius == true", textInput("radii", "Spheroid Radii (comma-separated):", "50,100")),
      selectInput("delivery_mode", "Delivery Mode (sa-mRNA and T cells):", choices=c("random","center","border"), selected="border"),
      numericInput("crowding_radius", "Crowding Radius (How far are neighbors):", 10, min=1, max=20),
      numericInput("crowding_beta", "Crowding Beta (Regulatory effect of neighbors):", 0.5, step=0.1),
      numericInput("ischemia_rate", "Ischemia Rate (Regulation by Ischemia):", 0.05, step=0.01),
      selectInput("ischemia_pattern", "Ischemia Pattern:", choices=c("center","none","periphery")),
      numericInput("tcell_kill_radius", "T Cell Kill Radius (How far T cells kill):", 20, min=1, max=50),
      numericInput("tcell_antigen_threshold", "T Cell Antigen Threshold (How much antigen required):", 2000, step=100),
      sliderInput("tcell_recognition_variability", "T Cell Recognition Variability:", min=0, max=1, value=0.1),
      numericInput("tcell_max_kills", "T Cell Max Kills per Hour (Time Step):", 0.2, min=0, max=10),
      sliderInput("infection_window", "Infection Window (h):", min=1, max=24, value=1),
      numericInput("n_replicates", "Number of Replicates:", 1, min=1, max=10, step=1),
      actionButton("run_sim", "Run Simulation"),
      hr(),
      downloadButton("downloadFig1", "Download Figure 1"),
      downloadButton("downloadFig2", "Download Figure 2"),
      downloadButton("downloadFig3", "Download Figure 3"),
      downloadButton("downloadFig4", "Download Figure 4")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Figure 1", plotOutput("fig1Plot", height="900px")),
        tabPanel("Figure 2", plotOutput("fig2Plot", height="900px")),
        tabPanel("Figure 3", plotOutput("fig3Plot", height="900px")),
        tabPanel("Figure 4", plotOutput("fig4Plot", height="900px"))
        
      )
    )
  )
)

# ===========================================================
# --- SERVER ---
# ===========================================================
server <- function(input, output, session) {
  
  sim_data <- eventReactive(input$run_sim, {
    set.seed(123)
    all_results <- list()
  
    # Parse user radii if enabled
    if (input$vary_radius) {
      radii <- as.numeric(unlist(strsplit(input$radii, ",")))
    } else {
      radii <- 50  # default single radius
    }
    
    for (replicate_id in 1:input$n_replicates) {
      for (spheroid_radius in radii) {
        
        n_cells <- input$n_cells
        spheroid_volume <- (4/3) * pi * spheroid_radius^3
        spheroid_density <- n_cells / spheroid_volume
        
        message(paste("Replicate", replicate_id,
                      "Radius", spheroid_radius,
                      "-> Volume =", round(spheroid_volume, 1),
                      "Density =", signif(spheroid_density, 3)))
        
        # Initialize positions and crowding
        positions <- initialize_positions(n_cells, geometry = input$geometry, L = spheroid_radius)
        cinfo <- compute_crowding(positions, 
                                  radius = input$crowding_radius, 
                                  beta = input$crowding_beta, 
                                  kernel = "gaussian", 
                                  mode = "inhibit")
        positions$crowding_factor <- cinfo$factor
        positions$neighbors <- cinfo$density
        
        # Initialize T cells
        tcell_positions <- initialize_tcells(input$n_tcells, geometry = input$geometry, L = spheroid_radius, placement = delivery_mode)
        
        # Run simulations
        reps <- lapply(1:nrow(positions), function(i){
          pos <- positions[i, ]
          simulate_cell(
            therapy = therapies$Replicon,
            dose = input$dose,
            cell_id = pos$cell_id,
            x = pos$x, y = pos$y, z = if("z" %in% names(pos)) pos$z else NULL,
            crowding_factor = pos$crowding_factor,
            ischemia_rate = input$ischemia_rate,
            ischemia_pattern = input$ischemia_pattern,
            spatial = TRUE,
            delivery_mode = input$delivery_mode,
            tcell_positions = tcell_positions,
            kill_radius = input$tcell_kill_radius,
            antigen_threshold = input$tcell_antigen_threshold,
            recognition_variability = input$tcell_recognition_variability,
            max_kills_per_tcell = input$tcell_max_kills,
            infection_window = input$infection_window,
            geometry = input$geometry
          )
        })
        
        df <- bind_rows(reps)
        df$dose <- input$dose
        df$replicate_id <- replicate_id
        df$spheroid_radius <- spheroid_radius
        df$spheroid_volume <- spheroid_volume
        df$spheroid_density <- spheroid_density
        df <- df |> left_join(positions[, c("cell_id","neighbors","crowding_factor")], by = "cell_id")
        
        all_results[[paste(replicate_id, spheroid_radius, sep = "_")]] <- df
      }
    }
    
    all_results <- bind_rows(all_results)

  })
  
  # --- Figure 1 ---
  fig1_obj <- reactive({
    df <- sim_data(); req(df)
    final_spatial <- df %>% filter(time==max(time)) %>% mutate(dist_center=sqrt(x^2 + y^2 + ifelse(!is.na(z), z^2,0)))
    agg_df <- df %>% group_by(time, neighbors) %>% summarise(mean_Antigen=mean(Antigen, na.rm=TRUE), .groups="drop")
    
    pA1 <- ggplot(final_spatial, aes(x=neighbors)) +
      geom_density(aes(group=replicate_id), color="grey70", alpha=0.5) +
      stat_density(geom="line", color="black", linewidth=1.2) +
      labs(title="Neighbors per Cell", x="Neighbor Count", y="Frequency")
    
    pB1 <- ggplot(agg_df, aes(x=time, y=mean_Antigen, group=neighbors, color=neighbors)) +
      geom_line(linewidth=1) + scale_color_viridis_c(option="plasma") +
      labs(title="GFP kinetics", x="Time (h)", y="GFP MFI")
    
    df_plot <- df %>% filter(time %in% c(0,6,12,18,24))
    pC1 <- ggplot(df_plot, aes(x=x, y=y, color=Antigen)) +
      geom_point(size=1.5) + scale_color_viridis(option="C") +
      facet_wrap(~time,ncol=5) + coord_fixed() + theme(panel.spacing=unit(2,"lines")) +
      labs(title="GFP spatial over time", x="X", y="Y")
    
    pD1 <- ggplot(final_spatial, aes(x=neighbors, y=Antigen)) +
      geom_point(alpha=0.3, color="green") + geom_smooth(method="loess", se=FALSE, color="black") +
      labs(title="GFP vs neighbors", x="Neighbor Count", y="GFP")
    
    pE1 <- ggplot(final_spatial, aes(x=dist_center, y=Antigen)) +
      geom_point(alpha=0.3, color="green") + geom_smooth(method="loess", se=FALSE, color="black") +
      labs(title="GFP vs distance", x="Distance", y="GFP")
    
    (pA1 | pB1) / (pC1) / (pD1 | pE1) + plot_annotation(tag_levels="A")
  })
  
  output$fig1Plot <- renderPlot({ fig1_obj() })
  
  # --- Figure 2 ---
  fig2_obj <- reactive({
    df <- sim_data(); req(df)
    df_24 <- df %>% filter(time==max(time))
    df_24$killed <- ifelse(df_24$alive,"Alive","Killed")
    
    pA2 <- ggplot(df_24, aes(x=x, y=y, color=killed)) +
      geom_point(size=2) + scale_color_manual(values=c("Alive"="darkgreen","Killed"="red")) +
      coord_fixed() + labs(title="Tumor Cell Status", x="X", y="Y")
    
    tcell_activity_filtered <- df_24 %>% filter(!is.na(killer_id)) %>%
      mutate(neighbors_int = round(neighbors)) %>%
      group_by(killer_id, neighbors_int) %>%
      summarise(kills = n(), .groups="drop")
    
    pB2 <- ggplot(tcell_activity_filtered, aes(x=factor(neighbors_int), y=kills)) +
      geom_jitter(width=0.2, alpha=0.1, color="darkblue") +
      geom_boxplot(fill="skyblue", alpha=0.5) +
      labs(title="T Cell Activity", x="Neighbor Count", y="Kills")+
      theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
    
    kill_24h_reps <- df_24 %>% mutate(neighbors_int = round(neighbors)) %>%
      group_by(neighbors_int) %>% summarise(frac_alive = mean(alive), .groups="drop")
    
    pC2 <- ggplot(kill_24h_reps, aes(x=neighbors_int, y=1-frac_alive)) +
      geom_point(size=2, alpha=0.3, position=position_jitter(width=0.2)) +
      geom_smooth() +
      labs(title="Tumor Cytotoxicity", x="Neighbor Count", y="Fraction Killed")
    
    (pA2 / pB2 / pC2) +
      plot_layout(heights = c(1, 1.5, 1)) +
      plot_annotation(tag_levels = "A")
  })
  
  output$fig2Plot <- renderPlot({ fig2_obj() })

  # --- Figure 3 ---
  fig3_obj <- reactive({
    df <- sim_data(); req(df)
    if (!"spheroid_radius" %in% names(df)) {
      showNotification("Figure 3: No spheroid radius data detected — enable 'Vary Radius' in settings.", type="warning")
      return(NULL)
    }
    
    final_spatial <- df %>%
      filter(time == max(time)) %>%
      mutate(dist_center = sqrt(x^2 + y^2 + ifelse(!is.na(z), z^2, 0)))
    
    agg_df <- df %>%
      group_by(time, neighbors, spheroid_radius) %>%
      summarise(mean_Antigen = mean(Antigen, na.rm = TRUE), .groups = "drop")
    
    # A: Neighbor distribution by radius
    pA3 <- ggplot(final_spatial, aes(x = neighbors, color = factor(spheroid_radius), fill = factor(spheroid_radius))) +
      geom_density(alpha = 0.3) +
      labs(title = "Neighbors per Cell by Spheroid Radius",
           x = "Neighbor Count", y = "Frequency", color = "Radius", fill = "Radius") +
      theme_minimal(base_size = 15)
    
    # B: GFP kinetics colored by radius
    df_summary <- df %>%
      group_by(time, spheroid_radius) %>%
      summarise(mean_Antigen = mean(Antigen, na.rm = TRUE),
                sd_Antigen = sd(Antigen, na.rm = TRUE),
                .groups = "drop")

    pB3 <- ggplot(df_summary, aes(x = time, y = mean_Antigen,
                                  color = factor(spheroid_radius),
                                  group = spheroid_radius)) +
      geom_line(size = 1.2) +
      geom_ribbon(aes(ymin = mean_Antigen - sd_Antigen,
                      ymax = mean_Antigen + sd_Antigen,
                      fill = factor(spheroid_radius)),
                  alpha = 0.5, color = NA) +
      labs(title = "Antigen Expression Kinetics by Spheroid Radius",
           x = "Time (h)",
           y = "Mean Antigen (GFP)",
           color = "Radius",
           fill = "Radius") +
      theme_minimal(base_size = 16)
    
    # C: Spatial plots over time, faceted by radius
    df_plot <- df %>% filter(time %in% c(0,6,12,18,24))
    pC3 <- ggplot(df_plot, aes(x = x, y = y, color = Antigen)) +
      geom_point(size = 1.2) +
      scale_color_viridis(option = "C") +
      facet_grid(spheroid_radius ~ time, labeller = label_both) +
      coord_fixed() +
      labs(title = "GFP Spatial Distribution over Time & Radius", x = "X", y = "Y") +
      theme_minimal(base_size = 14) +
      theme(panel.spacing = unit(1, "lines"))
    
    # D: GFP vs neighbors, colored by radius
    pD3 <- ggplot(final_spatial, aes(x = neighbors, y = Antigen,
                                     color = factor(spheroid_radius))) +
      geom_point(alpha = 0.5) +
      geom_smooth(se = FALSE) +
      labs(title = "GFP vs Neighbors", x = "Neighbor Count",
           y = "GFP", color = "Radius") +
      theme_minimal(base_size = 15)
    
    # E: GFP vs distance from center, colored by radius
    pE3 <- ggplot(final_spatial, aes(x = dist_center, y = Antigen,
                                     color = factor(spheroid_radius))) +
      geom_point(alpha = 0.5) +
      geom_smooth(se = FALSE) +
      labs(title = "GFP vs Distance", x = "Distance from Center",
           y = "GFP", color = "Radius") +
      theme_minimal(base_size = 15)
    
    (pA3 | pB3) / pC3 / (pD3 | pE3) + 
      plot_layout(heights = c(1, 1.5, 1)) +
      plot_annotation(tag_levels = "A")
  })
  
  output$fig3Plot <- renderPlot({ fig3_obj() })
  
  # --- Figure 4 ---
  fig4_obj <- reactive({
    df <- sim_data(); req(df)
    if (!"spheroid_radius" %in% names(df)) {
      showNotification("Figure 4: No spheroid radius data detected — enable 'Vary Radius' in settings.", type="warning")
      return(NULL)
    }
    
    df_24 <- df %>% filter(time == max(time))
    df_24$killed <- ifelse(df_24$alive, "Alive", "Killed")
    
    # A: Tumor cell status by radius
    pA4 <- ggplot(df_24, aes(x = x, y = y, color = killed)) +
      geom_point(size = 1.8, alpha = 0.7) +
      scale_color_manual(values = c("Alive" = "darkgreen", "Killed" = "red")) +
      coord_fixed() +
      facet_wrap(~ spheroid_radius, labeller = label_both) +
      labs(title = "Tumor Cell Status by Spheroid Radius", x = "X", y = "Y", color = "Status") +
      theme_minimal(base_size = 15)
    
    # B: T-cell activity, colored/faceted by radius
    tcell_activity <- df_24 %>%
      filter(!is.na(killer_id)) %>%
      mutate(neighbors_int = round(neighbors)) %>%
      group_by(killer_id, neighbors_int, spheroid_radius) %>%
      summarise(kills = n(), .groups = "drop")
    
    pB4 <- ggplot(tcell_activity, aes(x = factor(neighbors_int), y = kills,
                                      color = factor(spheroid_radius))) +
      geom_jitter(width = 0.2, alpha = 0.4) +
      geom_boxplot(alpha = 0.4, outlier.shape = NA) +
      theme(axis.text.x = element_text(angle = 90, hjust = 0.5))+
      labs(title = "T-Cell Activity by Spheroid Radius",
           x = "Neighbor Count", y = "Kills", color = "Radius") 
      
    
    # C: Cytotoxicity curves by radius
    kill_frac <- df_24 %>%
      mutate(neighbors_int = round(neighbors)) %>%
      group_by(neighbors_int, spheroid_radius) %>%
      summarise(frac_alive = mean(alive), .groups = "drop")
    
    pC4 <- ggplot(kill_frac, aes(x = neighbors_int, y = 1 - frac_alive,
                                 color = factor(spheroid_radius),
                                 group = spheroid_radius)) +
      geom_point(size = 2, alpha = 0.8, position = position_jitter(width = 0.2)) +
      geom_smooth(se = FALSE, linewidth = 1.2) +
      labs(title = "Tumor Cytotoxicity by Radius",
           x = "Neighbor Count", y = "Fraction Killed", color = "Radius") +
      theme_minimal(base_size = 15)
    
    (pA4 / pB4 / pC4) +
      plot_layout(heights = c(1, 1.5, 1)) +
      plot_annotation(tag_levels = "A")
  })
  
  output$fig4Plot <- renderPlot({ fig4_obj() })
  
  
  # --- Download handlers ---
  output$downloadFig1 <- downloadHandler(
    filename = function(){ "Figure1.png" },
    content = function(file){ ggsave(file, fig1_obj(), width=12, height=10, dpi=300) }
  )
  
  output$downloadFig2 <- downloadHandler(
    filename = function(){ "Figure2.png" },
    content = function(file){ ggsave(file, fig2_obj(), width=12, height=10, dpi=300) }
  )
  
  output$downloadFig3 <- downloadHandler(
    filename = function(){ "Figure3.png" },
    content = function(file){ ggsave(file, fig3_obj(), width=10, height=10, dpi=300) }
  )
  
  output$downloadFig4 <- downloadHandler(
    filename = function(){ "Figure4.png" },
    content = function(file){ ggsave(file, fig4_obj(), width=10, height=10, dpi=300) }
  )
  
}

# ===========================================================
# --- Run App ---
# ===========================================================
shinyApp(ui, server)
