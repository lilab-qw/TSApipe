kmerGeneration：对于癌症和正常样本，按顺序从原始样本生成 k-mer 数据库。
kmerFiltering：它的作用是对样本的kmer数据库进行过滤，过滤条件是：在control所用的mTEC生成的kmer文件中不存在，即count是0；并为每个k-mer添加一个ID，生成一个包含过滤后数据的新文件。这里涉及到一个LCOUNT，LCOUNT参数用于指定Jellyfish计算Cancer k-mer的最小出现次数，即只有在Cancer k-mer中出现次数不小于LCOUNT的k-mer才会被认为是表达的，当LCOUNT的值越高时，只有那些在样本中高度表达的k-mers才会被纳入考虑，因此最终生成的Cancer-specific k-mers数量会减少。相反，当LCOUNT的值较低时，更多的k-mers会被认为是表达的，这将导致生成更多的Cancer-specific k-mers。
一般而言，当数据质量较好时，可以使用较高的LCOUNT值，以生成更可靠和准确的Cancer-specific k-mers。反之，在数据质量较差的情况下，需要使用较低的LCOUNT值，以避免生成大量的噪声和误差。
kmerAssembly：