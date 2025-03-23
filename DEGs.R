library(DESeq2)
library(dplyr)

# 读取表达矩阵
mycounts <- read.csv("C:/Users/Lamarck/Desktop/mycounts.csv",row.names = 1)

# 查看读入的表达矩阵
head(mycounts)
dim(mycounts)

# 去除表达量全为0的行
mycounts_1 <- mycounts[rowSums(mycounts) != 0,]
#查看去除后的剩余的基因数目
dim(mycounts_1)

# 读取分组文件
mymeta <- read.csv("C:/Users/Lamarck/Desktop/mymeta.csv",stringsAsFactors = T)
mymeta
colnames(mycounts_1) == mymeta$id

dds <- DESeqDataSetFromMatrix(countData = mycounts_1, colData = mymeta,design = ~dex)
dds <- DESeq(dds)
res <- results(dds)

head(res)
res_1 <- data.frame(res)
head(res_1)

res_1 %>%
  mutate(group = case_when(
    log2FoldChange >= 1 & pvalue <= 0.05 ~ "UP",
    log2FoldChange <= -1 & pvalue <= 0.05 ~ "DOWN",
    TRUE ~ "NOT_CHANGE"
  )) -> res_2

table(res_2$group)


write.csv(res_2, file = "C:/Users/Lamarck/Desktop/diff_expr_result.csv", quote = F)
