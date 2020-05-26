calc_chi_u <- function(y, pair_ids, us){
  chis <- us
  for(i in 1:length(us)){
    xy <- list()
    for(k in 1:nrow(pair_ids)){
      xy[[k]] <- y[,pair_ids[k,]]
    }
    xy <- do.call(rbind, xy)
    xy <- xy[complete.cases(xy),]
    chis[i] <- taildep(xy[,1], xy[,2],u = us[i], type = "chi")
  }
  return(chis)
}

calc_chi_u_spatial <- function(y, coords, bins, us, bin_centers = NULL){
  D <- calc_dist(coords,coords)
  bl <- list()
  for (i in 1:(length(bins) - 1)) {
    bl[[i]] <- which(matrix(between(c(D), bins[i], bins[i +1]), 
                            nrow = nrow(D), ncol = ncol(D)), arr.ind = T)
  }
  if(is.null(bin_centers)){
    hs <- bins[-length(bins)] + diff(bins)/2
  }
  else{
    hs <- bin_centers
  }
  chi_list <- list()
  for(i in 1:length(bl)){
    if(nrow(bl[[i]]>0)){
      chi_list[[i]] <- data.frame(u = us, 
                                  chi = calc_chi_u(y, bl[[i]], us),
                                  bin = i,
                                  h = hs[i])
    }
  }
  do.call(rbind.data.frame, chi_list)
}