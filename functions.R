
#' Title
#'
#' @param edev.change 
#' @param tdev.change 
#' @param temp.change 
#' @param price.change 
#' @param pop.change 
#' @param dnetwork.change 
#' @param wqual.change 
#' @param edev.elasticity 
#' @param price.elasticity 
#' @param temp.elasticity 
#' @param dom.peak.factor 
#' @param ind.peak.factor 
#' @param year.sim 
#' @param year.ref 
#' @param nodes.data 
#' @param pipes.data 
#' @param global.output 
#'
#' @return
#' @export
#'
#' @examples
  simulateWDN <- function(edev.change = 0.0005,
                            tdev.change = 0.01,
                            temp.change = 0.04,
                            price.change = 0.01,
                            pop.change  = 0.01,
                            dnetwork.change = 0.01,
                            wqual.change = 0.005,
                            edev.elasticity = 1.0,
                            price.elasticity = -0.2,
                            temp.elasticity = 0.03,
                            dom.peak.factor = 1,
                            ind.peak.factor = 1,
                            year.sim = NULL,
                            year.ref = 2021,
                            nodes.data = NULL,
                            pipes.data = NULL,
                            global.output = TRUE)

  {

        # Demand & supply node indices
        dnodes <- nodes.data$id[which(nodes.data$type == "demand")]
        snodes <- nodes.data$id[which(nodes.data$type == "supply")]
        dnodes_num <- length(dnodes)
        snodes_num <- length(snodes)

        # Current year index
        yind <- year.sim - year.ref

        #:::::::::::::: DEMAND MODULE  :::::::::::::::::::::::::::::::::::::::::::::::

        # population growth rate per node (in the future, provide this as input)
        pop.change_pernode <- rep(pop.change, dnodes_num)

        # Demand adjustment factors
        edev_factor <- (1 + edev.change) ^ yind  #economic development factor for future year
        temp_factor <-  1 + temp.change * temp.elasticity * yind  # temperature factor for future year
        tdev_factor <- (1 + tdev.change) ^ yind #technological development factor for future year (for industrial demand)
        price_factor <- (1 + price.change * price.elasticity) ^ yind # price price factor for future year
        pop_factor <- (1 + pop.change_pernode) ^ yind #population factor for future year (per node)
        reuse_factor <- ifelse(yind < 23, 1, 0.90) # water reuse factor for future year
        dnetwork_factor <- (1 + dnetwork.change) ^ yind  # change in network demand percentage for future year

        # Redidential demand adjustment per node
        res_demand_adj <- temp_factor * price_factor * edev_factor * reuse_factor * pop_factor * dom.peak.factor

        # Industrial demand adjustment per node
        ind_demand_adj <- edev_factor * tdev_factor * temp_factor * dnetwork_factor * ind.peak.factor

        # Set demand per node
        demand_per_node_res <- nodes.data$discharge[dnodes] * (1-nodes.data$ind_ratio[dnodes]) * res_demand_adj
        demand_per_node_ind <- nodes.data$discharge[dnodes] * nodes.data$ind_ratio[dnodes] * ind_demand_adj

        # Total demand per node
        node_demand_update <- round(demand_per_node_res + demand_per_node_ind,2)
        nodes.data$discharge[dnodes] <- node_demand_update

        #:::::::::::::: SUPPLY MODULE  :::::::::::::::::::::::::::::::::::::::::::::::

        quality_factor <- (1 +  wqual.change) ^ yind
        node_supply_update = nodes.data$discharge[snodes]/quality_factor
        nodes.data$discharge[snodes] <- node_supply_update

        #::::::::::::: NETWORK MODULE ::::::::::::::::::::::::::::::::::::::::::::::::

        # simulate network
        res <- calculate_network(nodes.data, pipes.data, 1.0, 1)

        nodes_res = nodes.data %>%  as_tibble() %>%
            mutate(discharge = round(discharge, 2),
                   sim = abs(round(as.numeric(res[[2]]), 2)),
                   rel = round(100 * sim/discharge,0))

        pipes_res = pipes.data %>% as_tibble() %>%
            mutate(q_max = round(q_max, 2),
                   sim = abs(round(as.numeric(res[[3]]), 2)),
                   usage = round(100 * sim/q_max,0))


        if(global.output == TRUE) {

            # global_results
            result <- tibble(Parameter = c("T. Demand [m3]", "Available Sup. [m3]", "Allocated Sup. [m3]",
                                           "Reliability [%]", "Link usage [%]"),
                             value = NA)

            result$value[1] <- (nodes_res %>% filter(type == "demand") %>%
                                    filter(disuse == 0) %>% pull(discharge) %>% sum())

            result$value[2] <- (nodes_res %>% filter(type == "supply") %>%
                                    filter(disuse == 0) %>% pull(discharge) %>% sum())

            result$value[3] <- (nodes_res %>% filter(type == "supply") %>%
                                    filter(disuse == 0) %>% pull(sim) %>% sum())

            result$value[4] <- (nodes_res %>% filter(type == "demand") %>%
                                    filter(disuse == 0) %>%
                                    summarize(val = sum(sim)/sum(discharge) * 100) %>% pull(val))

            result$value[5] <- (pipes_res %>% filter(disuse == 0) %>%
                                    summarize(val = sum(sim)/sum(q_max) * 100) %>% pull(val))

            result$value <- as.numeric(round(result$value),2)

        } else {
            result <- NA
        }

        return(
            list(nodes = nodes_res, pipes = pipes_res, summary = result)
        )

    }

  
################################################################################
################################################################################
################################################################################
################################################################################
  
  
  
  
#' Title
#'
#' @param nodes.data 
#' @param pipes.data 
#' @param node.fill.var 
#' @param node.size.var 
#' @param edge.color.var 
#' @param edge.color.threshold 
#' @param edge.size.var 
#' @param background.map 
#'
#' @return
#' @export
#'
#' @examples
  visualizeWDN <- function(nodes.data = NULL,
                           pipes.data = NULL,
                           node.fill.var = NULL,
                           node.size.var = NULL,
                           edge.color.var = NULL,
                           edge.color.threshold = 95,
                           edge.size.var = NULL,
                           background.map = NULL)
  {

    require(dplyr)
    require(ggplot2)
    require(ggmap)
    require(ggnetwork)
    require(ggnewscale)
    require(scales)

    # Set ggplot2 theme
      ggtheme_network <- theme_light() +
          theme(legend.background = element_rect(fill = "white"),
                legend.key=element_blank(),
                legend.key.size = unit(0.6, "cm"))

      # Supply node aesthetics
      snode_color <- "blue"
      snode_fill <- "#c6dbef"
      snode_size <- 8
      snode_shape <- 22
      snode_stroke <- 1.2

      # Demand node aesthetics
      dnode_color <- "gray50"
      dnode_fill <- "gray50"
      dnode_size <- 8
      dnode_shape <- 21
      dnode_stroke <- 1


      # Edge/link aeasthetics
      edge_size <- 1
      edge_color <- c("0" = "red", "1" = "black", "2" = "green")

      # Axis labels
      xlab <- expression("Lon "(degree))
      ylab <- expression("Lat "(degree))

      # Prepare node/pipe tables
      nodesDF <- nodes.data
      pipesDF <- pipes.data %>%
          left_join(y = (nodesDF %>% select(start = id, lon, lat)), by = "start") %>%
          left_join(y = (nodesDF %>% select(end = id, lon_end = lon, lat_end = lat)), by = "end") %>%
          mutate(lon_mid = (lon + lon_end)/2,
                 lat_mid = (lat + lat_end)/2)


      if(is.null(node.fill.var)) {
          nodesDF$nodefill <- dnode_fill
      } else {
          nodesDF$nodefill <- nodesDF[[node.fill.var]]
      }

      if(is.null(node.size.var)) {
          nodesDF$nodesize <- dnode_size
      } else {
          nodesDF$nodesize <- nodesDF[[node.size.var]]
      }

      if(is.null(edge.color.var)) {
          pipesDF$edgecol <- 0
      } else {
          pipesDF$edgecol <- ifelse(pipesDF[[edge.color.var]] > edge.color.threshold, "0", "1")
      }

      pipesDF$edgecol <- ifelse(pipesDF[["disuse"]] > 0, "2", pipesDF$edgecol)

      if(is.null(edge.size.var)) {
          pipesDF$edgesize <- edge_size
      } else {
          pipesDF$edgesize <- pipesDF[[edge.size.var]]
      }

      edge_size_breaks <- c(100,300,500,700) #pipesDF %>% pull(diameter) %>% pretty(3)
      dnode_size_breaks <- nodesDF %>% filter(type == "demand") %>% pull(discharge) %>% pretty(5)

      nodesDF$nodefill[1] <- 100
      ##############################################################################

      # Main plot elements
      p <- ggmap(background.map, extend = "device")  +
          ggtheme_network +

          # Draw pipe network
          geom_edges(
              aes_string(x = "lon", y = "lat", xend = "lon_end", yend = "lat_end", color = "edgecol", size = "edgesize"),
              data = pipesDF) +
          scale_size(limits = range(edge_size_breaks), breaks = edge_size_breaks, range = c(0.6, 4),
                     guide = guide_legend(order = 3, title = "Pipe Diameter [mm]")) +
          scale_color_manual(values = edge_color, labels = c(" < 95", " >= 95", "NotInUse"),
                             guide = guide_legend(order = 4, title = "Pipe Use [%]",
                                                  override.aes = list(size = 3))) +
          new_scale("size") +

          # Draw supply nodes
          geom_nodes(
              mapping = aes_string(x = "lon", y = "lat"),
              data = filter(nodesDF, type == "supply"),
              fill = snode_fill, size = snode_size, stroke = snode_stroke,
              shape = snode_shape, color = snode_color) +

          # Draw demand nodes
          geom_nodes(
              aes_string(x = "lon", y = "lat", fill = "nodefill", size = "nodesize"),
              data = filter(nodesDF, type == "demand"),
              stroke = dnode_stroke, shape = dnode_shape, color = dnode_color) +

          scale_fill_steps2(low = "red", mid = "white", high = "white", midpoint = 99,
                            limits = c(0,100),
                            breaks = c(0,20, 40, 60, 80, 99, 100),
                            labels = c(0, 20, 40, 60, 80,"", 100),
                            guide = guide_colorbar(order = 1, title = "Reliability [%]",
                                                   frame.color = "black", barheight = unit(4, 'cm'))) +

          scale_size(limits = range(dnode_size_breaks), breaks = dnode_size_breaks, range = c(5,12),
                     guide = guide_legend(order = 2, title = "Demand [m3]")) +

          # Node labeling
          geom_text(aes(x = lon, y = lat, label = label), data = nodesDF, size=3.5) +

          # Set axis/guide labels
          labs(x = NULL, y = NULL)

      if(nrow(nodesDF %>% filter(disuse == 1)) > 0) {

          p <- p +

              geom_label(aes(x = lon, y = lat, label = "X"),
                         size = 5, alpha = 0.5, color = "white", fill = "green",
                         data = nodesDF %>% filter(disuse == 1))
      }


      if(nrow(pipesDF %>% filter(disuse == 1)) > 0) {

          p <- p +

              geom_label(aes(x = lon_mid, y = lat_mid, label = "X"),
                         size = 3, alpha = 0.5, color = "white", fill = "green",
                         data = pipesDF %>% filter(disuse == 1))
      }


      return(p)

  }
  