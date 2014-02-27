rankSVMcompare
==============

Support vector machines for ranking and comparing.

```s
install.packages(c("devtools","kernlab","quadprog","lpSolveAPI",
	           "ggplot2","directlabels"))
install.packages("quadmod", repos="http://r-forge.r-project.org")
library(devtools)
install_github("rankSVMcompare", "tdhock")
library(rankSVMcompare)
example(softCompareQP)
```
