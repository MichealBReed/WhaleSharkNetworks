############################################
#Whale Shark Analysis

#Author: Micheal B. Reed

###################################################################################
###Library
library('igraph')
library('SpiecEasi')
library('parallel')
library("plotly")
library('philentropy')


###################################################################################
#Functions

eigen_table <- function(igraph_object){
  sorted_eigen_names <- names(sort(eigen_centrality(igraph_object)$vector, decreasing = TRUE))
  sorted_eigen_values <- as.vector(sort(eigen_centrality(igraph_object)$vector, decreasing = TRUE))
  eigen_matrix <- matrix(data = sorted_eigen_values,
                         nrow = length(sorted_eigen_names),
                         ncol = 1)
  rownames(eigen_matrix) <- sorted_eigen_names
  return(eigen_matrix)
}


degree_table <- function(igraph_object){
  sorted_degree_names <- names(sort(degree(igraph_object), decreasing = TRUE))
  sorted_degree_values <- as.vector(sort(degree(igraph_object), decreasing = TRUE))
  degree_matrix <- matrix(data = sorted_degree_values,
                          nrow = length(sorted_degree_names),
                          ncol = 1)
  rownames(degree_matrix) <- sorted_degree_names
  return(degree_matrix)
}


between_table <- function(igraph_object){
  sorted_between_names <- names(sort(betweenness(igraph_object, directed = FALSE), decreasing = TRUE))
  sorted_between_values <- as.vector(sort(betweenness(igraph_object, directed = FALSE), decreasing = TRUE))
  between_matrix <- matrix(data = sorted_between_values,
                           nrow = length(sorted_between_names),
                           ncol = 1)
  rownames(between_matrix) <- sorted_between_names
  return(between_matrix)
}




taxa_cluster_table <- function(igraph_object){
  data_table <- data.frame(taxa = V(igraph_object)$name,
                           cluster = V(igraph_object)$membership)
  return(data_table)
}



preprocess_taxa <- function(data_file){
  "Function to make preprocessing more efficient. Assumes the first column are the Family names
  and the resulting columns are the taxa or funtions"
  phylo_names <- as.character(data_file[,1])
  
  inter_table <- data_file[,2:ncol(data_file)]
  
  rownames(inter_table) <- phylo_names
  
  if (length(grep('Avg', colnames(inter_table), ignore.case = TRUE)) > 0){
    
  average_column_idx <- grep('Avg', colnames(inter_table), ignore.case = TRUE)
  
  inter_table_2 <- inter_table[,colnames(inter_table)[-average_column_idx]]
  
  return(inter_table_2)
  
  }else{
    
    return(inter_table)
  }
  
}


preprocess_function <- function(data_file){
  
  function_names <- as.character(data_file[,1])
  
  inter_table <- data_file[,2:ncol(data_file)]
  
  rownames(inter_table) <- function_names
  
  inter_table_filtered <- inter_table[inter_table[ncol(inter_table)] > 100,]
  
  
  inter_table_filtered <- inter_table_filtered[,1:ncol(inter_table_filtered) - 1]
  
  return(inter_table_filtered)
  
}




normality_test <- function(igraph_object){
  return(shapiro.test(as.numeric(degree(igraph_object))))
}


qqplot <- function(igraph_object){
  degree_vector <- as.numeric(degree(igraph_object))
  qqnorm(degree_vector, pch = 1, frame = FALSE)
  qqline(degree_vector, col = 'violet', lwd = 2)
}

  


mcl_cluster <- function(igraph_object){
  mcl_result <- mcl(as_adjacency_matrix(igraph_object), addLoops = TRUE, inflation = 1.5, expansion = 2, max.iter = 200)
  return(mcl_result)
}



degree_viz <- function(igraph_object){
  random_model <- erdos.renyi.game(n=vcount(igraph_object),
                                   ecount(igraph_object),
                                   type = "gnm")
  
  barabasi_model <- barabasi.game(vcount(igraph_object),
                                  directed = FALSE)
  
  actual <- data.frame(degree = as.numeric(degree(igraph_object)))
  actual$model <- 'data'
  
  random <- data.frame(degree = as.numeric(degree(random_model)))
  random$model <- 'random'
  
  barabasi <- data.frame(degree = as.numeric(degree(barabasi_model)))
  barabasi$model <- 'scale-free'
  
  degrees <- rbind(actual,random,barabasi)
  
  p <- ggplot(degrees, aes(degree, fill = model)) + geom_density(alpha = 0.2)
  
  
  return(p)
  
}

####################################################################################################################
#Presence-Absence Analysis

taxa_pa_global_path <- "C:\\Users\\nemsr\\Desktop\\SDSU_projects\\WhaleSharkPaper\\Post_Sept_11\\Family_Data_9_11\\pres_abs_raw\\"

taxa_pa_file_list <- list(paste(taxa_pa_global_path, 'Cancun_family_ra.csv', sep = ''),
                          paste(taxa_pa_global_path, 'Ningaloo_family_ra.csv', sep = ''),
                          paste(taxa_pa_global_path, 'Philippines_family_ra.csv', sep = ''),
                          paste(taxa_pa_global_path, 'LaPaz_family_ra.csv', sep = ''))

taxa_raw_pa_data_frames <- lapply(taxa_pa_file_list, read.csv, header=TRUE)

taxa_processed_pa_data_frames <- lapply(taxa_raw_pa_data_frames, preprocess_taxa)

taxa_pa_binary_data_frames <- list(ifelse(taxa_processed_pa_data_frames[[1]] > 1,1,0),
                                   ifelse(taxa_processed_pa_data_frames[[2]] > 1,1,0),
                                   ifelse(taxa_processed_pa_data_frames[[3]] > 1,1,0),
                                   ifelse(taxa_processed_pa_data_frames[[4]] > 1,1,0)
)


taxa_pa_binary_data_frames_scale <- list(scale(taxa_pa_binary_data_frames[[1]],
                                               center = FALSE,
                                               scale = colSums(taxa_pa_binary_data_frames[[1]])),
                                         scale(taxa_pa_binary_data_frames[[2]],
                                               center = FALSE,
                                               scale = colSums(taxa_pa_binary_data_frames[[2]])),
                                         scale(taxa_pa_binary_data_frames[[3]],
                                               center = FALSE,
                                               scale = colSums(taxa_pa_binary_data_frames[[3]])),
                                         scale(taxa_pa_binary_data_frames[[4]],
                                               center = FALSE,
                                               scale = colSums(taxa_pa_binary_data_frames[[4]])))

Cancun_binary_se <- spiec.easi(t(taxa_pa_binary_data_frames_scale[[1]]),
                               method = 'mb',
                               lambda.min.ratio=1e-2,
                               nlambda=50
)


Ningaloo_binary_se <- spiec.easi(t(taxa_pa_binary_data_frames_scale[[2]]),
                               method = 'mb',
                               lambda.min.ratio=1e-2,
                               nlambda=50
)


Philippines_binary_se <- spiec.easi(t(taxa_pa_binary_data_frames_scale[[3]]),
                               method = 'mb',
                               lambda.min.ratio=1e-2,
                               nlambda=50
)


LaPaz_binary_se <- spiec.easi(t(taxa_pa_binary_data_frames_scale[[4]]),
                               method = 'mb',
                               lambda.min.ratio=1e-2,
                               nlambda=50
)

binary_se <- list(Cancun_binary_se,
                  Ningaloo_binary_se,
                  Philippines_binary_se,
                  LaPaz_binary_se)

binary_stability <- lapply(binary_se, getStability)

binary_network_list <- lapply(binary_se, getRefit)

binary_network_list <- lapply(binary_network_list, adj2igraph)




  
###################################################################################################################
#Preprocess our data so unneeded columns are not in our data set

family_global_path <- "C:\\Users\\nemsr\\Desktop\\SDSU_projects\\WhaleSharkPaper\\Post_Sept_11\\Family_Data_9_11\\100\\"

file_list_100 <- list(paste(family_global_path,"Cancun_family_100.csv",sep = ''),
                      paste(family_global_path,"Ningaloo_family_100.csv",sep = ''),
                      paste(family_global_path,"Philippines_family_100.csv",sep = ''),
                      paste(family_global_path,"LaPaz_family_100.csv",sep = ''))

raw_data_frames_100 <- lapply(file_list_100, read.csv, header=TRUE)

processed_data_frames_100 <- lapply(raw_data_frames_100, preprocess_taxa)

processed_data_frames_100 <- lapply(processed_data_frames_100, t)






Cancun_se_100 <- spiec.easi(processed_data_frames_100[[1]],
                        method = 'mb',
                        nlambda = 150,
                        lambda.min.ratio = 5e-2
                        )


Ningaloo_se_100 <- spiec.easi(processed_data_frames_100[[2]],
                          method = 'mb',
                          nlambda = 150,
                          lambda.min.ratio = 5e-2
                        )


Philippines_se_100 <- spiec.easi(processed_data_frames_100[[3]],
                             method = 'mb',
                             nlambda = 150)

LaPaz_se_100 <- spiec.easi(processed_data_frames_100[[4]],
                       method = 'mb',
                       nlambda = 150,
                       lambda.min.ratio = 5e-2
                       )



se_list_100 <- list(Cancun_se_100,
                Ningaloo_se_100,
                Philippines_se_100,
                LaPaz_se_100)

#Transform spiec-easi output into igraph objects
network_list_100 <- lapply(se_list_100, getRefit)

stability_list_100 <- lapply(se_list_100, getStability)

igraph_list_100 <- lapply(network_list_100, adj2igraph)


network_summary_100 <- data.frame(row.names = list('Cancun','Ningaloo','Philippines','LaPaz'),
                                  vertex_count = unlist(lapply(igraph_list_100, vcount)),
                                  edge_count = unlist(lapply(igraph_list_100, ecount)),
                                  diameter = unlist(lapply(igraph_list_100, diameter, directed=FALSE,weights=NA)),
                                  transitivity = unlist(lapply(igraph_list_100, transitivity, type='undirected', weights=NA))
)
                                  

#

V(igraph_list_100[[1]])$name <- colnames(processed_data_frames_100[[1]])

V(igraph_list_100[[2]])$name <- colnames(processed_data_frames_100[[2]])

V(igraph_list_100[[3]])$name <- colnames(processed_data_frames_100[[3]])

V(igraph_list_100[[4]])$name <- colnames(processed_data_frames_100[[4]])






                                 
                                 
              

#Visualize Degree Distributions
plot(degree.distribution(igraph_list_100[[1]]), type = 'b', ylab = 'Proportion', xlab = 'Degree')

plot(degree.distribution(igraph_list_100[[2]]), type = 'b', ylab = 'Proportion', xlab = 'Degree')

plot(degree.distribution(igraph_list_100[[3]]), type = 'b', ylab = 'Proportion', xlab = 'Degree')

plot(degree.distribution(igraph_list_100[[4]]), type = 'b', ylab = 'Proportion', xlab = 'Degree')



#Check the normality of Degree Distributions

normality_check_vector_100 <- unlist(lapply(igraph_list_100, normality_test))

#All Degree Distributions deviate from a normal distribution


#Random_network_vector

set.seed(123)
random_network_vector <- list(erdos.renyi.game(n=vcount(igraph_list_100[[1]]),
                                               ecount(igraph_list_100[[1]]),
                                               type = "gnm"),
                              erdos.renyi.game(n=vcount(igraph_list_100[[2]]),
                                               ecount(igraph_list_100[[2]]),
                                               type = "gnm"),
                              erdos.renyi.game(n=vcount(igraph_list_100[[3]]),
                                               ecount(igraph_list_100[[3]]),
                                               type = "gnm"),
                              erdos.renyi.game(n=vcount(igraph_list_100[[4]]),
                                               ecount(igraph_list_100[[4]]),
                                               type = "gnm"))


random_degree_distribution_vector <- lapply(random_network_vector, degree.distribution)


barabasi_network_vector <- list(barabasi.game(vcount(igraph_list_100[[1]]), directed = FALSE),
                                barabasi.game(vcount(igraph_list_100[[2]]), directed = FALSE),
                                barabasi.game(vcount(igraph_list_100[[3]]), directed = FALSE),
                                barabasi.game(vcount(igraph_list_100[[4]]), directed = FALSE)
                                )

barabasi_degree_distribution_vector <- lapply(barabasi_network_vector, degree.distribution)

#Visualize degree distributions vs random network of similar dimensions
plot(degree.distribution(igraph_list_100[[1]]), type = 'b', ylab = 'Proportion', xlab = 'Degree', xlim = c(0,7), ylim = c(0,0.6))
points(0:(length(random_degree_distribution_vector[[1]])-1), random_degree_distribution_vector[[1]], col="red" , type='b')
points(0:(length(barabasi_degree_distribution_vector[[1]])-1), barabasi_degree_distribution_vector[[1]], col="blue" , type='b')


plot(degree.distribution(igraph_list_100[[2]]), type = 'b', ylab = 'Proportion',
     xlab = 'Degree', xlim = c(0,8), ylim = c(0, max(degree.distribution(igraph_list_100[[2]]))))
points(0:(length(random_degree_distribution_vector[[2]])-1), random_degree_distribution_vector[[2]], col="red" , type='b')


plot(degree.distribution(igraph_list_100[[3]]), type = 'b', ylab = 'Proportion', xlab = 'Degree', xlim = c(0,10))
points(0:(length(random_degree_distribution_vector[[3]])-1), random_degree_distribution_vector[[3]], col="red" , type='b')

plot(degree.distribution(igraph_list_100[[4]]), type = 'b', ylab = 'Proportion', xlab = 'Degree', xlim = c(0,8))
points(0:(length(random_degree_distribution_vector[[4]])-1), random_degree_distribution_vector[[4]], col="red" , type='b')



#Visuzalize the degree distributions from all locations in comparison with analogous random and scale-free networks
degree_viz(igraph_list_100[[1]])

degree_viz(igraph_list_100[[2]])

degree_viz(igraph_list_100[[3]])

degree_viz(igraph_list_100[[4]])



#######################################################################
#Markov Clustering

#MCL clustering

CL <- makeCluster(detectCores())



clusterEvalQ(CL,{
  library(SpiecEasi)
  library(MCL)
  library(igraph)
})



#I have 4 cpus on my laptop and so I use these 4 cores to run the mcl process on four of our datasets
mcl_easi_list <- parLapply(CL,
                           igraph_list_100,
                           mcl_cluster)

#close the cluster
stopCluster(CL)

V(igraph_list_100[[1]])$membership <- mcl_easi_list[[1]]$Cluster

V(igraph_list_100[[2]])$membership <- mcl_easi_list[[2]]$Cluster

V(igraph_list_100[[3]])$membership <- mcl_easi_list[[3]]$Cluster

V(igraph_list_100[[4]])$membership <- mcl_easi_list[[4]]$Cluster



#Summary table having location, taxa, betweenness, eigen_centrality and cluster membership

taxa_summary_table <- data.frame(location = c(rep('Cancun', vcount(igraph_list_100[[1]])),
                                              rep('Ningaloo', vcount(igraph_list_100[[2]])),
                                              rep('Philippines', vcount(igraph_list_100[[3]])),
                                              rep('LaPaz', vcount(igraph_list_100[[4]]))),
                                 taxa = c(V(igraph_list_100[[1]])$name,
                                          V(igraph_list_100[[2]])$name,
                                          V(igraph_list_100[[3]])$name,
                                          V(igraph_list_100[[4]])$name),
                                 
                                 between_cen = c(betweenness(igraph_list_100[[1]], directed = FALSE, weights = NA),
                                                 betweenness(igraph_list_100[[2]], directed = FALSE, weights = NA),
                                                 betweenness(igraph_list_100[[3]], directed = FALSE, weights = NA),
                                                 betweenness(igraph_list_100[[4]], directed = FALSE, weights = NA)),
                                 
                                 eigen_cen = c(as.numeric(eigen_centrality(igraph_list_100[[1]])$vector),
                                               as.numeric(eigen_centrality(igraph_list_100[[2]])$vector),
                                               as.numeric(eigen_centrality(igraph_list_100[[3]])$vector),
                                               as.numeric(eigen_centrality(igraph_list_100[[4]])$vector)),
                                 
                                 cluster_mem = c(paste('Ca.', V(igraph_list_100[[1]])$membership, sep = ''),
                                                 paste('Ni.', V(igraph_list_100[[2]])$membership, sep = ''),
                                                 paste('Ph.', V(igraph_list_100[[3]])$membership, sep = ''),
                                                 paste('LaP.', V(igraph_list_100[[4]])$membership, sep = ''))
)




#################################################################################################
#Functions Network Analysis

function_global_path <- "C:\\Users\\nemsr\\Desktop\\SDSU_projects\\WhaleSharkPaper\\Post_Sept_11\\Functional_Data_9_29\\"

Cancun_function_file_100 <- read.csv(paste(function_global_path,'Cancun_sept_29.csv', sep = ''), header = TRUE)

Ningaloo_function_file_100 <- read.csv(paste(function_global_path,'Ningaloo_sept_29.csv', sep = ''), header = TRUE)

Philippines_function_file_100 <- read.csv(paste(function_global_path,'Philippines_sept_29.csv', sep = ''), header = TRUE)

LaPaz_function_file_100 <- read.csv(paste(function_global_path,'LaPaz_sept_29.csv', sep = ''), header = TRUE)



raw_function_file_list <- list(Cancun_function_file_100,
                           Ningaloo_function_file_100,
                           Philippines_function_file_100,
                           LaPaz_function_file_100)


processed_function_files <- lapply(raw_function_file_list, preprocess_function)

processed_function_files <- lapply(processed_function_files, as.matrix)

processed_function_files <- lapply(processed_function_files, floor)

processed_function_files <- lapply(processed_function_files, t)

##############################################################################################
#Function Spiec-Easi

Cancun_function_se <- spiec.easi(processed_function_files[[1]],
                                                  method = 'mb',
                                                  nlambda = 150
                                 )
          


Ningaloo_function_se <- spiec.easi(processed_function_files[[2]],
                                   method = 'mb',
                                   nlambda = 200
                                   
)


Philippines_function_se <- spiec.easi(processed_function_files[[3]],
                                      method = 'mb',
                                      nlambda = 150
                                      
)


LaPaz_function_se <- spiec.easi(processed_function_files[[4]],
                                method = 'mb',
                                nlambda = 150
                                
)



function_se_list <- list(Cancun_function_se,
                         Ningaloo_function_se,
                         Philippines_function_se,
                         LaPaz_function_se)

function_stability_list_100 <- lapply(function_se_list, getStability)


function_network_list_100 <- lapply(function_se_list, getRefit)

function_stability_list_100 <- lapply(function_se_list, getStability)

function_igraph_list_100 <- lapply(function_network_list_100, adj2igraph)



V(function_igraph_list_100[[1]])$name <- colnames(processed_function_files[[1]])

V(function_igraph_list_100[[2]])$name <- colnames(processed_function_files[[2]])

V(function_igraph_list_100[[3]])$name <- colnames(processed_function_files[[3]])

V(function_igraph_list_100[[4]])$name <- colnames(processed_function_files[[4]])





#Markov Clustering 

#MCL clustering

CL <- makeCluster(detectCores())



clusterEvalQ(CL,{
  library(SpiecEasi)
  library(MCL)
  library(igraph)
})



#I have 4 cpus on my laptop and so I use these 4 cores to run the mcl process on four of our datasets
function_mcl_easi_list <- parLapply(CL,
                           function_igraph_list_100,
                           mcl_cluster)

#close the cluster
stopCluster(CL)



V(function_igraph_list_100[[1]])$membership <- function_mcl_easi_list[[1]]$Cluster

V(function_igraph_list_100[[2]])$membership <- function_mcl_easi_list[[2]]$Cluster

V(function_igraph_list_100[[3]])$membership <- function_mcl_easi_list[[3]]$Cluster

V(function_igraph_list_100[[4]])$membership <- function_mcl_easi_list[[4]]$Cluster



function_summary_table <- data.frame(location = c(rep('Cancun', vcount(function_igraph_list_100[[1]])),
                                              rep('Ningaloo', vcount(function_igraph_list_100[[2]])),
                                              rep('Philippines', vcount(function_igraph_list_100[[3]])),
                                              rep('LaPaz', vcount(function_igraph_list_100[[4]]))),
                                 taxa = c(V(function_igraph_list_100[[1]])$name,
                                          V(function_igraph_list_100[[2]])$name,
                                          V(function_igraph_list_100[[3]])$name,
                                          V(function_igraph_list_100[[4]])$name),
                                 
                                 between_cen = c(betweenness(function_igraph_list_100[[1]], directed = FALSE, weights = NA),
                                                 betweenness(function_igraph_list_100[[2]], directed = FALSE, weights = NA),
                                                 betweenness(function_igraph_list_100[[3]], directed = FALSE, weights = NA),
                                                 betweenness(function_igraph_list_100[[4]], directed = FALSE, weights = NA)),
                                 
                                 eigen_cen = c(as.numeric(eigen_centrality(function_igraph_list_100[[1]])$vector),
                                               as.numeric(eigen_centrality(function_igraph_list_100[[2]])$vector),
                                               as.numeric(eigen_centrality(function_igraph_list_100[[3]])$vector),
                                               as.numeric(eigen_centrality(function_igraph_list_100[[4]])$vector)),
                                 
                                 cluster_mem = c(paste('Ca.', V(function_igraph_list_100[[1]])$membership, sep = ''),
                                                 paste('Ni.', V(function_igraph_list_100[[2]])$membership, sep = ''),
                                                 paste('Ph.', V(function_igraph_list_100[[3]])$membership, sep = ''),
                                                 paste('LaP.', V(function_igraph_list_100[[4]])$membership, sep = ''))
                                 
                                
)


degree_viz(function_igraph_list_100[[1]])

degree_viz(igraph_list_100[[2]])

degree_viz(igraph_list_100[[3]])

degree_viz(igraph_list_100[[4]])


function_network_summary_100 <- data.frame(row.names = list('Cancun','Ningaloo','Philippines','LaPaz'),
                                  vertex_count = unlist(lapply(function_igraph_list_100, vcount)),
                                  edge_count = unlist(lapply(function_igraph_list_100, ecount)),
                                  diameter = unlist(lapply(function_igraph_list_100, diameter, directed=FALSE,weights=NA)),
                                  transitivity = unlist(lapply(function_igraph_list_100, transitivity, type='undirected', weights=NA))
)

stability_list_100 <- lapply(se_list_100, getStability)



##########################################################
#Graph Visualization

plot(igraph_list_100[[1]],
     vertex.color = V(igraph_list_100[[1]])$membership,
     vertex.label = NA,
     layout = layout.fruchterman.reingold(igraph_list_100[[1]]),
     cex = 5)

plot(igraph_list_100[[2]],
     vertex.color = V(igraph_list_100[[2]])$membership,
     vertex.label = NA,
     layout = layout.fruchterman.reingold(igraph_list_100[[2]]),
     cex = 5)


