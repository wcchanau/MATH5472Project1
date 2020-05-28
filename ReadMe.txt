VB_clustering.Rdata
variable	description
mu		true cluster center
K		number of cluster
n		size of training dataset
n_testing	size of testing dataset
new_order	used in re-arranging the cluster index
sigma2		variance of mu
c		cluster index of training dataset
c_testing	cluster index of testing dataset
y		observed value of training dataset
y_testing	observed value of testing dataset
phi_testing	the probability of belonging to each cluster based on the result from training
parameter	training result
--m		estimated posterior cluster center (after reorder)
--s2		estimated posterior covariance matrix of cluster center (after reorder)
--pi		the probability of belonging to each cluster 
--L		evidence lower bound

