## Script for getting nucleosome distance matrix ## 
## Execute the file like this: Rscript folder

source("/orozco/homes/pluto/jwalther/Programs/R/image_scale.R")

#Input single file

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).xyz", call.=FALSE)
		      }

#out_folder = "/home/MuG/MuG_Chromatin_sampling/analysis"
#file="/orozco/homes/pluto/jwalther/Programs/Master_Manuel/Manuel_Chromatin_equ_structure/src_test/output_nucl/nucl_000000.xyz"


folder = args[1]
#folder="/orozco/homes/jwalther/Programs/MuG/MuG_Chromatin_sampling/src_test/output_nucl"
files = list.files(folder, full.names=TRUE)

out_f = args[2]

out_folder=sprintf("%s/output", out_f)

loop_op = function(ind)
{
file = files[ind]	# Input file is file of x,y,z coordinates of nucleosomes

b = read.table(file, skip=2)
b = b[,-1]
d = as.matrix(dist(b, method = "euclidean", diag = FALSE, upper = TRUE))

mean_dist_df = as.data.frame(d)

#### Ensemble parameters ####

# Radius of gyration
names(b) = c("x","y","z")
N = length(b[,1])
r_mean = unlist(lapply(1:3, function(i) mean(as.numeric(b[,i]))))

vec = lapply(1:length(b[,1]), function(i) (b[i,] - r_mean)*(b[i,] - r_mean))
vec_df = do.call(rbind.data.frame, vec)
N1 = length(vec)
squ_norm_vec = unlist(lapply(1:length(vec), function(x) sum(vec_df[x,])))

r_g = sqrt(1/N1 * sum(squ_norm_vec))
# END Radius of gyration


# Fit a function to each dimension x,y,z (Perisic et al. (2010))
nu = c(1:N)

acc=c(10,10,10)
if(N<10) acc=c(N-1,N-1,N-1)
#acc=c(4,4,4)

model_x = lm(b$x ~ poly(nu,acc[1], raw=TRUE))
model_y = lm(b$y ~ poly(nu,acc[2], raw=TRUE))
model_z = lm(b$z ~ poly(nu,acc[3], raw=TRUE))

pred_x = predict(model_x, newdata = data.frame(x=b$x))
pred_y = predict(model_y, newdata = data.frame(y=b$y))
pred_z = predict(model_z, newdata = data.frame(z=b$z))

##### packing ratio
pred_d = cbind(pred_x,pred_y,pred_z)
d = as.matrix(dist(pred_d, method = "euclidean", diag = FALSE, upper = FALSE))
cont_len = sum(unlist(lapply(nu[-length(nu)], function(x) d[x,x+1])))	# Take every nucleosome
#cont_len1 = sum(unlist(lapply(seq(1,N-2,2), function(x) d[x,x+2])))	# Take every second nucleosome

pack_ratio = N/(cont_len/10/11)
#pack_ratio1 = N/(cont_len1/10/11)

##### end to end distance
dis = as.matrix(dist(b, method = "euclidean", diag = FALSE, upper = FALSE))
end_to_end_d = dis[length(b[,1]),1]


##### Mean fiber radius
co = model_x$coef
inter = seq(1,N,0.2)	# Adjust according to accuracy (0.2 is fine for long fibers)
predacc_x = colSums(do.call("rbind",lapply(1:length(co), function(i) co[i]*inter^(i-1))))
co = model_y$coef
predacc_y = colSums(do.call("rbind",lapply(1:length(co), function(i) co[i]*inter^(i-1))))
co = model_z$coef
predacc_z = colSums(do.call("rbind",lapply(1:length(co), function(i) co[i]*inter^(i-1))))


pred = unname(cbind(predacc_x,predacc_y,predacc_z))
min_dist = function(i)
{
buf = rbind(unname(as.numeric(b[i,])),unname(pred))
d = as.matrix(dist(buf, method = "euclidean", diag = FALSE, upper = FALSE))
min_d = min(d[2:length(d[,1]),1])
return(min_d)
}

all_min_d = unlist(lapply(1:N, min_dist))
mean_all_min = mean(all_min_d)/10 + 5.5	# converted from A to nm and added 5.5nm for radius of nucleosome (as in Perisic et al (2010))
d_all_min = 2*mean_all_min	# fiber diameter

#### END Ensemble parameters ####

## Output for VRE ##
result = round(c(cont_len,end_to_end_d,pack_ratio,d_all_min,r_g),2)
names(result) = c("Contour length (in nm)","End-to-end distance (in nm)","Packing ratio","Fiber diameter (in nm)","Radius of gyration (in nm)")

return(list(mean_dist_df,result))
}


res = lapply(1:length(files), loop_op)

mat = lapply(res, function(x) x[[1]])
mean_dist = Reduce("+", mat) / length(mat)


ens_pars = do.call("rbind", lapply(res, function(x) x[[2]]))
mean = round(colMeans(ens_pars),2)
sd = round(apply(ens_pars, 2, sd),2)
means_sds_df = rbind(mean,sd)
means_sds_df[,c(1,2,4,5)] = means_sds_df[,c(1,2,4,5)]/10

write.csv(means_sds_df, file=sprintf("%s/create_traj_out.csv",out_folder))

#d = as.matrix(dist(b, method = "euclidean", diag = FALSE, upper = TRUE))

d = as.matrix(mean_dist)

n=200
rc1 = heat.colors(n)
#rc1 = rainbow(30)
white = "#FFFFFFFF"
rc1 = c(rc1,rep(white,n/0.4))
e = as.numeric(d)

min=min(e[!(e == min(e))])
max=max(e)
subimg_w = c(10,2)
img_height = 1024
img_width = (subimg_w[1]+subimg_w[2])/subimg_w[1]*img_height
png(sprintf("%s/create_traj01.png",out_folder),width=img_width, height=img_height, bg="white")
layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(subimg_w[1],subimg_w[2]), heights=c(10,10))
#layout.show(4)
par(mar=c(1,1,1,1))
image(d, zlim=c(min,max), col=rc1, xaxt="n", yaxt="n")
par(mar=c(1,1,1,5))
image.scale(d, zlim=c(min,max), col=rc1,horiz=FALSE, xaxt="n", yaxt="n")
axis(4, cex.axis=2)
dev.off()



####### Neighboring distances of nucleosomes #######

mean_dist_df = as.data.frame(d)

neigh_dist = function(i)
{
x = as.numeric(mean_dist_df[i,])

get_dist = function(j)
{
if((i-j) > 0) a=x[i-j] else a=NA
if((i+j) < length(x)) b=x[i+j] else b=NA
c(a,b)
}
res = do.call("rbind",lapply(1:length(x),get_dist))
}

neigh_d = lapply(1:length(mean_dist_df[,1]),neigh_dist)

neigh_d_cmb = do.call("rbind",neigh_d)

group = function(i)
{
a = seq(i,length(neigh_d_cmb[,1]),length(neigh_d[[1]][,1]))
x = neigh_d_cmb[a,]
x = as.numeric(x)
x = x[!is.na(x)]
x = unique(x)
return(x)
}
grouped_neighd = lapply(1:length(neigh_d[[1]][,1]),group)

grouped_mean_sd = do.call("rbind",lapply(grouped_neighd, function(x) c(mean(x),sd(x))))


bord=length(grouped_mean_sd[,1])-1
x = c(1:bord)
mi = min(log(grouped_mean_sd[1:bord,1]))-0.1
ma = max(log(grouped_mean_sd[1:bord,1]))+0.1

png(sprintf("%s/create_traj02.png",out_folder),width=600, height=600, res=100, bg="white")
plot(x,log(grouped_mean_sd[1:bord,1]), ylim=c(mi,ma), main=sprintf("Decay of internucleosomal interactions with distance"), xlab="Distance (in # of nucleosomes)",ylab="log(distance (in A))",col="blue")
dev.off()



