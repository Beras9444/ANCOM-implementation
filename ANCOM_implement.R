Before_table=read.csv(<First time point>)
After_table=read.csv(<Second time point>)

#make a pseudocount
for(x in seq(2,ncol(After_table),1)){
  for(y in seq(1,3,1)){
    After_table[y,x]=After_table[y,x]+1
  }
}
for(x in seq(2,ncol(Before_table),1)){
  for(y in seq(1,3,1)){
    Before_table[y,x]=as.numeric(Before_table[y,x])+1
  }
}


#Adding a column at the end to contain sum of rows
colu=ncol(Before_table)

for (q in seq(1,3)){
  sumv1=c()
  for (p in seq(2,colu)){
    sumv1=c(sumv1,as.numeric(Before_table[q,p]))
  }
  Before_table[q,colu+1]=sum(sumv1)
}

colu=ncol(After_table)

for (q in seq(1,3)){
  sumv2=c()
  for (p in seq(2,colu)){
    sumv2=c(sumv2,as.numeric(After_table[q,p]))
  }
  After_table[q,colu+1]=sum(sumv2)
}

#convert to relative abundance
for(y in seq(2,colu)){
  for(w in seq(1,3)){
    Before_table[w,y]=((as.numeric(Before_table[w,y]))/(as.numeric(Before_table[w,colu+1])))*100
    After_table[w,y]=((as.numeric(After_table[w,y]))/(as.numeric(After_table[w,colu+1])))*100
  }
}

#Find log of mean Relative Abundance
Before_table[4,] = c('log of mean of i', as.numeric(apply(Before_table[,-1],2, function(x) log(mean(x)))))
After_table[4,] = c('log of mean of i', as.vector(apply(After_table[,-1],2, function(x) log(mean(x)))))


#Find a reference taxon
start=(as.numeric(Before_table[4,2])-as.numeric(After_table[4,2]))^2
ref=2
for(t in seq(2,625)){
  val<-((as.numeric(After_table[4,t]))-(as.numeric(Before_table[4,t])))^2
  if(val < start){
    ref<-t
    start<-val
    }
}



#Convert to log ratios over reference taxon
for(y in seq(2,colu)){
  for(w in seq(1,3)){
    if (y==ref){
      next
    }
    Before_table[w,y]=(log(as.numeric(Before_table[w,y])))-(log(as.numeric(Before_table[w,ref])))
    After_table[w,y]=log(as.numeric(After_table[w,y]))-log(as.numeric(After_table[w,ref]))
  }
}

#Find log mean ratios over mean reference in a new row
Before_table[5,] = c('log of mean of i/mean of ref', as.numeric(apply(Before_table[4,c(-1)],2, function(x) as.numeric(x)-as.numeric(Before_table[4,ref]))))
After_table[5,] = c('log of mean of i/mean of ref', as.vector(apply(After_table[4,c(-1)],2, function(x) as.numeric(x)-as.numeric(After_table[4,ref]))))
  



#Perform a t test for all taxons, and report the taxons with p values below 0.05
i_list=c()
p_list=c()
for (n in seq(2,colu)){
  arr1=c(as.numeric(Before_table[c(-4,-5),n]))
  arr2=c(as.numeric(After_table[c(-4,-5),n]))
  test_r=t.test(arr1,arr2)
  if(n==ref){
    i_list=c(i_list,"ref")
  }else{
    if(test_r$p.value<0.05){
    i_list=c(i_list,n)
    Before_table[6,n]=test_r$p.value
    }
  }

}

#Convert list to numeric form
i_list=as.numeric(i_list)



#Report list with significantly different taxons based on log ratios, followed by their difference in log ratios of means over reference

taxa_list2=c()
for(r in i_list){
  if(is.na(r)){
    next
  }
  taxa_list2=c(taxa_list2,colnames(Before_table)[r])
}

#Remove NA from reference taxon from list
list2=i_list[-which(is.na(i_list))]

#Find ratio of log of mean taxon ratios over ref taxon of each significantly differing taxon 
#Also making a list for p values
change=c()
p_value=c()
for (k in list2){
  change=c(change,(as.numeric(Before_table[5,k])-as.numeric(After_table[5,k])))
  p_value=c(p_value,as.numeric(Before_table[6,k]))
}


#Make list of all significantly differing taxons names using indices frim previous list
change_list<-colnames(Before_table)[list2]

#Make list of significantly differing taxons, followed by p value and their multiple of change based on log ratio of mean of taxon over mean of ref taxon
final_list=c()
for (r in seq(1,length(change_list))){
  final_list=c(final_list,c(change_list[r],p_value[r],change[r]))
}

write(final_list,file=<Write location>)
