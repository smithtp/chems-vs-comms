# take the 8 chemicals, randomly order them to produce 2 and 4 chem mixes that make up the 8 chem mixture.

library(dplyr)
library(tidyr)
library(ggnetwork)

stressor_vector <- c("A", "C", "D", "G", "I", "M", "O", "T")

stressor_sample <- sample(stressor_vector)

pairs <- list(paste(sort(stressor_sample[1:2]), collapse = ""),
              paste(sort(stressor_sample[3:4]), collapse = ""),
              paste(sort(stressor_sample[5:6]), collapse = ""),
              paste(sort(stressor_sample[7:8]), collapse = ""))

quads <- list(paste(sort(stressor_sample[1:4]), collapse = ""),
              paste(sort(stressor_sample[5:8]), collapse = ""))


# could do something cool like plot it on to the network (showing only the nodes and lines we've selected)
# or maybe show all the nodes, but only the edges we want, and have names only for the nodes we want.

# step one, build all nodes and edges as previously

chem_matrix <- read.csv("data/chemical-matrix-plus.csv")

chem_matrix <- chem_matrix %>%
  arrange(Complexity, ChemCode)


# create edges
# edges will be from the complexity level below, which contain all but one of the chems in the upper node

from_vector <- c()
to_vector <- c()

for(i in 2:8){
  from_options <- chem_matrix[chem_matrix$Complexity == i-1,]$ChemCode
  to_options <- chem_matrix[chem_matrix$Complexity == i,]$ChemCode
  for(j in 1:length(from_options)){
    for(k in 1:length(to_options)){
      condition <- unlist(strsplit(from_options[j], split = "")) %in% unlist(strsplit(to_options[k], split = ""))
      if(!(FALSE %in% condition)){
        from_vector <- c(from_vector, from_options[j])
        to_vector <- c(to_vector, to_options[k])
      }
    }
  }
  
}


all_edges <- data.frame(source = from_vector, destination = to_vector)

all_nodes <- data.frame(label = unique(c(all_edges$source, all_edges$destination)))
all_nodes$id <- rownames(all_nodes)

nodes <- data.frame(id = all_nodes$id, label = all_nodes$label)

# now need to link the edges to the ID labels in nodes
edges <- all_edges %>% 
  left_join(nodes, by = c("source" = "label")) %>% 
  rename(from = id)

edges <- edges %>% 
  left_join(nodes, by = c("destination" = "label")) %>% 
  rename(to = id)

edges <- select(edges, from, to)

# this lays out exactly where the nodes should be:
nodes$x <- c(seq(21,49, by = 4), seq(8,62, by = 2), seq(2,68, by = 1.2), seq(1,70), seq(2,68, by = 1.2), seq(8,62, by = 2), seq(21,49, by = 4), 35) 
nodes$y <- c(rep(1, 8), rep(2, 28), rep(3, 56), rep(4, 70), rep(5, 56), rep(6, 28), rep(7, 8), 8)

# create the network from the edges and nodes
chem_network_gg <- ggnetwork(edges, layout = as.matrix(nodes[,c("x", "y")]))

# add the chemical codes labels back onto the matrix dataframe
chem_network_gg <- left_join(chem_network_gg, nodes[,c("id", "label")], by = c("vertex.names" = "id"))

######################
# now the fancy part #
######################

# try #1 - rename node labels so only the ones we've selected are shown
sample_network <- chem_network_gg[chem_network_gg$label %in% c(pairs, quads, stressor_sample, "ACDGIMOT"),]

ggplot(chem_network_gg, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_nodelabel_repel(aes(label = label), fontface = "bold", box.padding = unit(1, "lines")) +
  #geom_edges(color = "grey50") +
  geom_nodes() +
  theme_blank()

ggplot(sample_network, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_nodelabel_repel(aes(label = label), fontface = "bold", box.padding = unit(1, "lines")) +
  geom_edges(color = "grey50") +
  geom_nodes() +
  theme_blank()


# try # 2 - keep all nodes, but only label/colour the ones we've selected

sample_network <- chem_network_gg %>%
  mutate(sample_label = NA)
sample_network[sample_network$label %in%  c(pairs, quads, stressor_sample, "ACDGIMOT"),]$sample_label <- 
  sample_network[sample_network$label %in%  c(pairs, quads, stressor_sample, "ACDGIMOT"),]$label

# draw some new edges corresponding to the sampled data
# the network has "vertex names" which seem to just be the nodes
# the x and y correspond to where the vertex sits
# the xend and yend correspond to a line from the node (x,y), finishing on the row above (xend,yend)
# would be simpler to just draw them touching the points they go to, i.e. xend/yend = x/y of the edge finish

# create edges from scratch
new_network <- sample_network[0,]

# function to create edges
for(i in 1:8){
  temp_network <- sample_network %>% filter(sample_label == stressor_sample[i]) %>% slice(1)
  # create new xend/yend based on where the first pair is
  temp_end <- sample_network %>% filter(sample_label == pairs[[ceiling(i/2)]]) %>% slice(1)
  # replace xend and yend in first df with x and y in 2nd df
  temp_network$xend <- temp_end$x
  temp_network$yend <- temp_end$y
  # bind into dataframe
  new_network <- rbind(new_network, temp_network)
}

for(i in 1:4){
  temp_network <- sample_network %>% filter(sample_label == pairs[[i]]) %>% slice(1)
  # create new xend/yend based on where the first pair is
  temp_end <- sample_network %>% filter(sample_label == quads[[ceiling(i/2)]]) %>% slice(1)
  # replace xend and yend in first df with x and y in 2nd df
  temp_network$xend <- temp_end$x
  temp_network$yend <- temp_end$y
  # bind into dataframe
  new_network <- rbind(new_network, temp_network)
}

for(i in 1:2){
  temp_network <- sample_network %>% filter(sample_label == quads[[i]]) %>% slice(1)
  # create new xend/yend based on where the first pair is
  temp_end <- sample_network %>% filter(sample_label == "ACDGIMOT") %>% slice(1)
  # replace xend and yend in first df with x and y in 2nd df
  temp_network$xend <- temp_end$x
  temp_network$yend <- temp_end$y
  # bind into dataframe
  new_network <- rbind(new_network, temp_network)
}

ggplot(sample_network, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(data = new_network, color = "red") +
  geom_nodes() +
  geom_nodelabel(aes(label = sample_label), fontface = "bold", box.padding = unit(1, "lines")) +
  # geom_nodelabel_repel(aes(label = sample_label), fontface = "bold", box.padding = unit(1, "lines")) +
  theme_blank()



# turn this all into a single function that we can run
sampling_function <- function(plot = TRUE){
  # first sample the stressors
  stressor_sample <- sample(stressor_vector)
  
  # pull out the pairs
  pairs <- list(paste(sort(stressor_sample[1:2]), collapse = ""),
                paste(sort(stressor_sample[3:4]), collapse = ""),
                paste(sort(stressor_sample[5:6]), collapse = ""),
                paste(sort(stressor_sample[7:8]), collapse = ""))
  
  # pull out the 4-combinations
  quads <- list(paste(sort(stressor_sample[1:4]), collapse = ""),
                paste(sort(stressor_sample[5:8]), collapse = ""))
  
  # alter the network so only the sampled parts are named
  sample_network <- chem_network_gg %>%
    mutate(sample_label = NA)
  sample_network[sample_network$label %in%  c(pairs, quads, stressor_sample, "ACDGIMOT"),]$sample_label <- 
    sample_network[sample_network$label %in%  c(pairs, quads, stressor_sample, "ACDGIMOT"),]$label
  
  # create edges from scratch
  new_network <- sample_network[0,]
  
  # function to create edges
  for(i in 1:8){
    temp_network <- sample_network %>% filter(sample_label == stressor_sample[i]) %>% slice(1)
    # create new xend/yend based on where the first pair is
    temp_end <- sample_network %>% filter(sample_label == pairs[[ceiling(i/2)]]) %>% slice(1)
    # replace xend and yend in first df with x and y in 2nd df
    temp_network$xend <- temp_end$x
    temp_network$yend <- temp_end$y
    # bind into dataframe
    new_network <- rbind(new_network, temp_network)
  }
  
  for(i in 1:4){
    temp_network <- sample_network %>% filter(sample_label == pairs[[i]]) %>% slice(1)
    # create new xend/yend based on where the first pair is
    temp_end <- sample_network %>% filter(sample_label == quads[[ceiling(i/2)]]) %>% slice(1)
    # replace xend and yend in first df with x and y in 2nd df
    temp_network$xend <- temp_end$x
    temp_network$yend <- temp_end$y
    # bind into dataframe
    new_network <- rbind(new_network, temp_network)
  }
  
  for(i in 1:2){
    temp_network <- sample_network %>% filter(sample_label == quads[[i]]) %>% slice(1)
    # create new xend/yend based on where the first pair is
    temp_end <- sample_network %>% filter(sample_label == "ACDGIMOT") %>% slice(1)
    # replace xend and yend in first df with x and y in 2nd df
    temp_network$xend <- temp_end$x
    temp_network$yend <- temp_end$y
    # bind into dataframe
    new_network <- rbind(new_network, temp_network)
  }
  
  p <- ggplot(sample_network, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(data = new_network, color = "red") +
    geom_nodes() +
    geom_nodelabel(aes(label = sample_label), fontface = "bold") +
    # geom_nodelabel_repel(aes(label = sample_label), fontface = "bold", box.padding = unit(1, "lines")) +
    theme_blank()
  
  if(plot == TRUE){
    print(p)
  }
  
  return(c(pairs, quads))
}

sampling_function()


###############################################################################
# write a function to make 2 --unique-- pyramids from each community experiment
# and save the output in a useful way
# also randomizing which 2 of the 4 replicates we sample for each pyramid


# list of communities
communities_list <- c("R_ID_7.16_1_R10(1)_TSS",
                      "R_ID_7.16_4_R10(1)_TSS",
                      "R_ID_7.16_5_R10(1)_TSS",
                      "R_ID_7.16_6_R5(1)_TSS",
                      "R_ID_7.16_7_R1(1)_TSS",
                      "R_ID_7.16_8_R5(1)_TSS",
                      "R_ID_7.16_9_R1(1)_TSS",
                      "R_ID_7.16_10_R10(1)_TSS",
                      "R_ID_7.16_11_R11(1)_TSS",
                      "R_ID_7.16_12_R5(1)_TSS)")

communities_numbers <- seq(1, 10)

# we need to sample all the single wells and the 8-way mixture regardless of anything else
static_rows <- chem_matrix %>% filter(Complexity %in% c(1,8))

# place to put the results
community_sampling_df <- data.frame()

for(i in seq_along(communities_numbers)){
  # We need to sample the starting community
  starting_community <- data.frame(Community.number = communities_numbers[i],
                                   Community.name = communities_list[i],
                                   ChemCode = "Frozen Starting Community")
  
  community_sampling_df <- bind_rows(community_sampling_df, starting_community)
  
  # randomly sample 2 control wells
  all_controls <- bind_rows(chem_matrix %>% filter(Complexity == 0) %>%
                              mutate(Plate.sample = 5,
                                     Well.sample = paste(Plate.sample, Well, sep = " ")),
                            chem_matrix %>% filter(Complexity == 0) %>%
                              mutate(Plate.sample = 10,
                                     Well.sample = paste(Plate.sample, Well, sep = " ")),
                            chem_matrix %>% filter(Complexity == 0) %>%
                              mutate(Plate.sample = 15,
                                     Well.sample = paste(Plate.sample, Well, sep = " ")),
                            chem_matrix %>% filter(Complexity == 0) %>%
                              mutate(Plate.sample = 20,
                                     Well.sample = paste(Plate.sample, Well, sep = " "))) %>%
    select(-X)
  
  control_sample <- all_controls[sample(nrow(all_controls), 2),] %>%
    mutate(Community.number = communities_numbers[i],
           Community.name = communities_list[i])
  
  community_sampling_df <- bind_rows(community_sampling_df, control_sample)
  
  # then loop through the static rows and randomly sample 2 of the 4 reps
  for(j in 1:nrow(static_rows)){
    static_row <- static_rows[j,]
    plate_num <- static_row$Plate
    plate_reps <- c(0, 5, 10, 15)
    
    plate_vector <- plate_num + plate_reps
    plate_sample <- sample(plate_vector, size = 2)
    
    for(k in 1:2){
      selected_row <- static_row %>%
        select(-X) %>%
        mutate(Community.number = communities_numbers[i],
               Community.name = communities_list[i],
               Plate.sample = plate_sample[k],
               Well.sample = paste(Plate.sample, Well, sep = " "),
               Pyramid = 1)
      
      community_sampling_df <- bind_rows(community_sampling_df, selected_row)
      
    }
  }
  # now we've done the static data for this community, do the sampling of the pyramids
  # sample the first pyramid
  set_1 <- sampling_function(plot = FALSE)
  # pull out the dataframe rows
  set_1_rows <- chem_matrix %>% filter(ChemCode %in% set_1)
  # loop through the rows and randomly sample 2 of the 4 reps
  for(p in 1:nrow(set_1_rows)){
    sampled_row <- set_1_rows[p,]
    plate_num <- sampled_row$Plate
    plate_reps <- c(0, 5, 10, 15)
    
    plate_vector <- plate_num + plate_reps
    plate_sample <- sample(plate_vector, size = 2)
    
    for(k in 1:2){
      selected_row <- sampled_row %>%
        select(-X) %>%
        mutate(Community.number = communities_numbers[i],
               Community.name = communities_list[i],
               Plate.sample = plate_sample[k],
               Well.sample = paste(Plate.sample, Well, sep = " "),
               Pyramid = 1)
      
      community_sampling_df <- bind_rows(community_sampling_df, selected_row)
      
    }
  }
  
  # sample second pyramid
  # put this inside a logical thing to check the two pyramids are unique
  repeat {
    # do something
    set_2 <- sampling_function(plot = FALSE)
    set_2_rows <- chem_matrix %>% filter(ChemCode %in% set_2)
    # exit if the condition is met
    if (!any(set_1 %in% set_2)) break
  }
  
  # loop through the rows and randomly sample 2 of the 4 reps
  for(q in 1:nrow(set_2_rows)){
    sampled_row <- set_2_rows[q,]
    plate_num <- sampled_row$Plate
    plate_reps <- c(0, 5, 10, 15)
    
    plate_vector <- plate_num + plate_reps
    plate_sample <- sample(plate_vector, size = 2)
    
    for(k in 1:2){
      selected_row <- sampled_row %>%
        select(-X) %>%
        mutate(Community.number = communities_numbers[i],
               Community.name = communities_list[i],
               Plate.sample = plate_sample[k],
               Well.sample = paste(Plate.sample, Well, sep = " "),
               Pyramid = 2)
      
      community_sampling_df <- bind_rows(community_sampling_df, selected_row)
      
    }
  }
}

head(community_sampling_df)

# wonderful, now write this out

write.csv(community_sampling_df, "communities-to-extract.csv", row.names = FALSE)
