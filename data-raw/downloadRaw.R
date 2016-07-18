
## Download raw data
temp <- tempfile()
download.file("http://www3.nd.edu/~jli9/mseq/data_top100.zip",temp)
filenames = unzip(temp, list = TRUE)$Name
filenames = filenames[grep('.csv', filenames)]
lnames = substr(filenames,13,14)
for(jj in 1:length(filenames)) assign(lnames[jj], read.csv(unz(temp,filenames[jj])))
unlink(temp)

devtools::use_data(b1,b2,b3,g1,g2,w1,w2,w3,internal=TRUE)

## Process these files
l = 20
r = 19
b1 = mseq::expData(b1, l, r)[,-1]
b2 = mseq::expData(b2, l, r)[,-1]
b3 = mseq::expData(b3, l, r)[,-1]
g1 = mseq::expData(g1, l, r)[,-1]
g2 = mseq::expData(g2, l, r)[,-1]
w1 = mseq::expData(w1, l, r)[,-1]
w2 = mseq::expData(w2, l, r)[,-1]
w3 = mseq::expData(w3, l, r)[,-1]

devtools::use_data(b1,b2,b3,g1,g2,w1,w2,w3)
