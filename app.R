  library(shiny)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(shinythemes)
  library(gganimate)
  library(gifski)
  library(png)
  
  
  # ============================================================
  #   Simulación SIRD → SIRVD con vacuna imperfecta (modelo nuevo)
  # ============================================================
  
  # ------------------------------------------------------------
  # Simulación del modelo determinista. 
  simulate_sird_sirdv <- function(
    
    # Parámetros iniciales
    N = 15000,
    I0 = 10,
    R0 = 0,
    D0 = 0,
    V0 = 0,
    days = 150,
    dt = 0.1,
    
    beta0 = 0.6,
    gamma = 1/10,
    mu = 0.005,
    
    day_confinamiento = 30,
    confinamiento_factor = 0.6,
    
    day_vacuna = 70,
    nu = 0.02,            # tasa de vacunación (probabilidad por susceptible por día)
    eficacia = 0.8,       # epsilon
    theta = 1             # velocidad resolución vacuna
  ){
    
    # tiempo
    t <- seq(0, days, by = dt)
    
    # variables
    S <- rep(0, length(t))
    V <- rep(0, length(t))
    I <- rep(0, length(t))
    R <- rep(0, length(t))
    D <- rep(0, length(t))
    
    # condiciones iniciales
    S[1] <- N - I0 - R0 - D0 - V0
    I[1] <- I0
    R[1] <- R0
    D[1] <- D0
    V[1] <- V0
    
    # Beta dependiente del confinamiento
    beta_t <- function(day){
      if(day < day_confinamiento) return(beta0)
      else return(beta0 * confinamiento_factor)
    }
    
    # Tasa de vacunación dependiente del día
    nu_t <- function(day){
      if(day < day_vacuna) return(0)
      else return(nu)
    }
    
    # ===============================
    #      Simulación paso a paso
    # ===============================
    for(k in 1:(length(t)-1)){
      day <- t[k]
      
      beta <- beta_t(day)
      nu_day <- nu_t(day)
      
      N_actual <- S[k] + V[k] + I[k] + R[k]
      
      lambda <- beta * I[k] / max(N_actual, 1e-9)
      
      # Ecuaciones DEL MODELO NUEVO
      dS <- -lambda * S[k] - nu_day * S[k] + (1 - eficacia) * theta * V[k]
      dV <- nu_day * S[k] - theta * V[k]
      dI <- lambda * S[k] - gamma * I[k] - mu * I[k]
      dR <- gamma * I[k] + eficacia * theta * V[k]
      dD <- mu * I[k]
      
      # Integración Euler
      S[k+1] <- S[k] + dS * dt
      V[k+1] <- V[k] + dV * dt
      I[k+1] <- I[k] + dI * dt
      R[k+1] <- R[k] + dR * dt
      D[k+1] <- D[k] + dD * dt
    }
    
    # salida
    data.frame(
      dia = t,
      S = S,
      V = V,
      I = I,
      R = R,
      D = D
    )
  }
  
  # ------------------------------------------------------------
  # Simulación para el modelo estocástico
  
  simulate_sird_sirdv_stochastic <- function(
    N = 15000,
    I0 = 10,
    R0 = 0,
    D0 = 0,
    V0 = 0,
    days = 150,
    dt = 0.1,
    
    beta0 = 0.6,
    gamma = 1/10,
    mu = 0.005,
    
    day_confinamiento = 30,
    confinamiento_factor = 0.6,
    
    day_vacuna = 70,
    nu = 0.02,
    eficacia = 0.8,
    theta = 1
  ){
    
    t <- seq(0, days, by = dt)
    n <- length(t)
    
    S <- rep(0, n)
    V <- rep(0, n)
    I <- rep(0, n)
    R <- rep(0, n)
    D <- rep(0, n)
    
    S[1] <- N - I0 - R0 - D0 - V0
    I[1] <- I0
    R[1] <- R0
    D[1] <- D0
    V[1] <- V0
    
    # Beta dependiente del tiempo
    beta_t <- function(day){
      if(day < day_confinamiento) beta0 else beta0 * confinamiento_factor
    }
    
    # Tasa de vacunación dependiente del tiempo
    nu_t <- function(day){
      if(day < day_vacuna) 0 else nu
    }
    
    for(k in 1:(n-1)){
      day <- t[k]
      
      beta <- beta_t(day)
      nu_day <- nu_t(day)
      
      N_actual <- max(S[k] + V[k] + I[k] + R[k], 1e-9)
      lambda <- beta * I[k] / N_actual
      
      # --- Convert to transition probabilities ---
      p_inf  <- min(lambda * dt, 1)
      p_vacc <- min(nu_day * dt, 1)
      p_wane <- min((1 - eficacia) * theta * dt, 1)
      p_full <- min(eficacia * theta * dt, 1)
      p_rec  <- min(gamma * dt, 1)
      p_die  <- min(mu * dt, 1)
      
      # --- Stochastic transitions ---
      
      # S transitions
      new_inf  <- rbinom(1, S[k], p_inf)
      S_after_inf <- S[k] - new_inf
      
      new_vacc <- rbinom(1, S_after_inf, p_vacc)
      
      # V transitions
      waned <- rbinom(1, V[k], p_wane)
      V_after_wane <- V[k] - waned
      
      vaccinated_recovered <- rbinom(1, V_after_wane, p_full)
      
      # I transitions
      recov <- rbinom(1, I[k], p_rec)
      death <- rbinom(1, I[k] - recov, p_die)
      
      # --- Update compartments ---
      S[k+1] <- S[k] - new_inf - new_vacc + waned
      V[k+1] <- V[k] + new_vacc - waned - vaccinated_recovered
      I[k+1] <- I[k] + new_inf - recov - death
      R[k+1] <- R[k] + recov + vaccinated_recovered
      D[k+1] <- D[k] + death
    }
    
    data.frame(
      dia = t,
      S = S,
      V = V,
      I = I,
      R = R,
      D = D
    )
  }
  # ------------------------------------------------------------
  # Simulación para estocástico con múltiples líneas
  # ------------------------------------------------------------
  
  simulate_multiple_stochastic <- function(reps = 20, ...) {
    sims <- lapply(1:reps, function(i) {
      df <- simulate_sird_sirdv_stochastic(...)
      df$sim_id <- i
      df
    })
    dplyr::bind_rows(sims)
  }
  
  # ------------------------------------------------------------
  # Función de ayuda para los estilos
  # ------------------------------------------------------------
  
  theme_cyber <- function() {
    theme_minimal(base_size = 15) +
      theme(
        plot.background = element_rect(fill = "#0d0d0d", color = NA),
        panel.background = element_rect(fill = "#0d0d0d", color = NA),
        legend.background = element_rect(fill = "#0d0d0d", color = NA),
        panel.grid.major = element_line(color = "#222222"),
        panel.grid.minor = element_line(color = "#111111"),
        text = element_text(color = "#cccccc"),
        axis.text = element_text(color = "#cccccc"),
        axis.title = element_text(color = "#eeeeee"),
        plot.title = element_text(color = "#00ffea", face = "bold")
      )
  }
  
  neon_palette <- c(
    "S" = "#00eaff",
    "I" = "#ffee00",
    "R" = "#00ff88",
    "V" = "#ff00ff",
    "D" = "#ff3355",
    "mean" = "#ffffff"
  )
  
  # ------------------------------------------------------------
  # Función para la animación
  # ------------------------------------------------------------
  animate_curve <- function(df, title = "Animación", filename){
    
    df_long <- df %>% pivot_longer(cols = -dia)
    
    p <- ggplot(df_long, aes(dia, value, color = name)) +
      geom_line(size = 1.4) +
      scale_color_manual(values = neon_palette) +
      labs(title = title, x = "Día", y = "Población") +
      theme_cyber() +
      transition_reveal(dia)
    
    anim <- animate(
      p,
      nframes = 200,
      fps = 20,
      width = 600,
      height = 400,
      renderer = gifski_renderer(loop = TRUE)
    )
    
    anim_save(filename, animation = anim)
    filename
  }
  
  
  # ============================================================
  #                          UI
  # ============================================================
  
  ui <- fluidPage(
    theme = shinytheme("cyborg"),
    titlePanel("Simulación Epidemiológica SIRD → SIRDV Vibrio Silvestris"),
    
    sidebarLayout(
      sidebarPanel(
        
        numericInput("N", "Población total (N):", 50000),
        numericInput("I0", "Infectados iniciales (I0):", 3),
        numericInput("days", "Duración (días):", 80),
        
        hr(),
        sliderInput("beta0", "β (transmisión):", 0.1, 1, 0.4, step=0.01),
        sliderInput("gamma", "γ (recuperación):", 0.05, 0.3, 0.08, step=0.01),
        sliderInput("mu", "μ (mortalidad):", 0, 0.02, 0.005, step=0.0005),
        
        hr(),
        numericInput("day_confin", "Día confinamiento:", 30),
        sliderInput("conf_factor", "Reducción β:", 0.1, 1, 0.35),
        
        hr(),
        numericInput("day_vac", "Día inicio vacuna:", 40),
        sliderInput("nu", "Tasa vacunación (ν):", 0, 0.2, 0.2, step=0.005),
        sliderInput("efi", "Eficacia vacuna (ε):", 0, 1, 0.5),
        sliderInput("theta", "Velocidad efecto vacuna (θ):", 0.1, 3, 1, step=0.1),
        
        hr(),
        checkboxInput("toggle_confin", "Activar confinamiento", TRUE),
        checkboxInput("toggle_vacuna", "Activar vacuna", TRUE),
        
      ),
      
      mainPanel(
        tabsetPanel(
          tabPanel("Determinístico",
                   actionButton("run", "Ejecutar modelo"),
                   plotOutput("plot_det"),
                   tableOutput("table_det"),
                   actionButton("anim_det_btn", "Generar Animación Determinística"),
                   imageOutput("anim_det")
          ),
          tabPanel("Estocástico",
                   actionButton("run_stoch", "Ejecutar modelo estocástico"),
                   plotOutput("plot_stoch"),
                   tableOutput("table_stoch"),
                   hr(),
                   h3("Simulación Multitrayectoria (Estocástico)"),
                   numericInput("replicas", "Número de trayectorias:", 20, min = 2, max = 200),
                   actionButton("run_stoch_multi", "Ejecutar múltiples trayectorias"),
                   plotOutput("plot_stoch_multi", height = "600px"),
                   tableOutput("table_stoch"),
                   actionButton("anim_stoch_btn", "Generar Animación Estocástica"),
                   imageOutput("anim_stoch")
                   
          ))))
  )
  
  # ============================================================
  #                        SERVER
  # ============================================================
  
  server <- function(input, output){
    
    # Modelo determinista
    data_det <- eventReactive(input$run, {
      
      # Read inputs
      day_confin <- input$day_confin
      conf_factor <- input$conf_factor
      day_vac <- input$day_vac
      nu_val <- input$nu
      
      # ---- Apply toggles ----
      if (!input$toggle_confin) {
        day_confin <- Inf
        conf_factor <- 1
      }
      
      if (!input$toggle_vacuna) {
        nu_val <- 0
        day_vac <- Inf
      }
      
      simulate_sird_sirdv(
        N = input$N,
        I0 = input$I0,
        days = input$days,
        beta0 = input$beta0,
        gamma = input$gamma,
        mu = input$mu,
        day_confinamiento = day_confin,
        confinamiento_factor = conf_factor,
        day_vacuna = day_vac,
        nu = nu_val,
        eficacia = input$efi,
        theta = input$theta
      )
    })
    
    
    # Modelo estocástico
    data_stoch <- eventReactive(input$run_stoch, {
      
      day_confin <- input$day_confin
      conf_factor <- input$conf_factor
      day_vac <- input$day_vac
      nu_val <- input$nu
      
      if (!input$toggle_confin) {
        day_confin <- Inf
        conf_factor <- 1
      }
      
      if (!input$toggle_vacuna) {
        nu_val <- 0
        day_vac <- Inf
      }
      
      simulate_sird_sirdv_stochastic(
        N = input$N,
        I0 = input$I0,
        days = input$days,
        beta0 = input$beta0,
        gamma = input$gamma,
        mu = input$mu,
        day_confinamiento = day_confin,
        confinamiento_factor = conf_factor,
        day_vacuna = day_vac,
        nu = nu_val,
        eficacia = input$efi,
        theta = input$theta
      )
    })
    
    
    ##################################
    # PLOTS
    ##################################
    # Plots determinista
    output$plot_det <- renderPlot({
      df <- data_det()
      df_long <- df %>% pivot_longer(cols = -dia)
      
      ggplot(df_long, aes(dia, value, color = name)) +
        geom_line(size = 1.3) +
        geom_vline(xintercept = input$day_confin,
                   linetype="dashed", color="#ff3355", size=1) +
        geom_vline(xintercept = input$day_vac,
                   linetype="dotted", color="#ff00ff", size=1) +
        scale_color_manual(values = neon_palette) +
        theme_cyber() +
        labs(
          title = "Modelo Determinístico SIRDV",
          x = "Día", y = "Población"
        )
    })
    
    # Tabla modelo determinista
    output$table_det <- renderTable({
      head(data_det(), 10)
    })
    
    # Animación modelo determinístico
    observeEvent(input$anim_det_btn, {
      output$anim_det <- renderImage({
        df <- data_det()
        outfile <- tempfile(fileext = ".gif")
        
        animate_curve(
          df,
          title = "Animación: Modelo Determinístico",
          filename = outfile
        )
        
        list(src = outfile, contentType = "image/gif")
      }, deleteFile = TRUE)
    })
    
    
    # Plots estocástico
    output$plot_stoch <- renderPlot({
      df <- data_stoch()
      df_long <- df %>% pivot_longer(cols = -dia)
      
      ggplot(df_long, aes(dia, value, color = name)) +
        geom_line(size = 1.3) +
        geom_vline(xintercept = input$day_confin,
                   linetype="dashed", color="#ff3355", size=1) +
        geom_vline(xintercept = input$day_vac,
                   linetype="dotted", color="#ff00ff", size=1) +
        scale_color_manual(values = neon_palette) +
        theme_cyber() +
        labs(
          title = "Modelo Estocástico SIRDV (1 Simulación)",
          x = "Día", y = "Población"
        )
    })
  
    # Tabla modelo estocástico
    output$table_stoch <- renderTable({
      head(data_stoch(), 10)
    })
    
    # Animación modelo estocástico
    observeEvent(input$anim_stoch_btn, {
      output$anim_stoch <- renderImage({
        df <- data_stoch()
        outfile <- tempfile(fileext = ".gif")
        
        animate_curve(
          df,
          title = "Animación: Modelo Estocástico",
          filename = outfile
        )
        
        list(src = outfile, contentType = "image/gif")
      }, deleteFile = TRUE)
    })
    
    # --------------------------------------------------------------------------
    # Modelo estocástico Multilinea
    # --------------------------------------------------------------------------
    observeEvent(input$run_stoch_multi, {
      output$plot_stoch_multi <- renderPlot({
        
        reps <- input$replicas
        
        # Run multiple trajectories
        df_all <- simulate_multiple_stochastic(
          reps = reps,
          N = input$N,
          I0 = input$I0,
          days = input$days,
          beta0 = input$beta0,
          gamma = input$gamma,
          mu = input$mu,
          day_confinamiento = input$day_confin,
          confinamiento_factor = input$conf_factor,
          day_vacuna = input$day_vac,
          nu = input$nu,
          eficacia = input$efi,
          theta = input$theta
        )
        
        # Compute mean across simulations
        df_mean <- df_all %>%
          group_by(dia) %>%
          summarise(across(c(S, I, R, V, D), mean)) %>%
          mutate(sim_id = "mean")
        
        # Reshape all data
        df_all_long <- df_all %>%
          pivot_longer(cols = c(S,I,R,V,D), names_to = "comp", values_to = "value")
        
        df_mean_long <- df_mean %>%
          pivot_longer(cols = c(S,I,R,V,D), names_to = "comp", values_to = "value")
        
        # Plot
        ggplot() +
          # All stochastic trajectories
          geom_line(data = df_all_long,
                    aes(x = dia, y = value, group = interaction(sim_id, comp),
                        color = comp),
                    alpha = 0.25, size = 0.6) +
          # Mean curve
          geom_line(data = df_mean_long,
                    aes(x = dia, y = value, color = comp),
                    size = 1.7) +
          scale_color_manual(values = neon_palette) +
          labs(
            title = paste("Trayectorias Estocásticas (", reps, " simulaciones)", sep = ""),
            x = "Día",
            y = "Población"
          ) +
          theme_cyber()
      })
    })
  }
  
  # ============================================================
  #                    Ejecutar aplicación
  # ============================================================
  
  shinyApp(ui, server)
