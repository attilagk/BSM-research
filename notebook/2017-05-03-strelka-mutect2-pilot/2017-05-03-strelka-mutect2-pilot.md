

## Scaling with segment length


```r
runt <- read.csv("results/runtimes.csv")
```


```r
#tpar <- trellis.par.get()
#trellis.par.set("superpose.line", list(lty = 2))
arg <- list(xlab = "length (bases)", ylab = "real time",
            labels = c("1s", "7.8s", "1m", "7.8m", "1h", "7.8h", "2.5d", "2.8w", "5m"),
            at.lin = 0:8 / 2, at.log = 60 ^ c(0:8 / 2))
p1 <- xyplot(runtime.strelka + runtime.mutect2 ~ seglen, data = runt, type = "b", xlab = arg$xlab, ylab = arg$ylab)
p2 <- xyplot(log(runtime.strelka, base = 60) + log(runtime.mutect2, base = 60) ~ log10(seglen), data = runt, type = c("p", "r"), xlab = arg$xlab, ylab = arg$ylab, scales = list(y = list(at = arg$at.lin, labels = arg$labels)), grid = TRUE)
print(p1, split = c(1, 1, 2, 1), more = TRUE)
print(p2, split = c(2, 1, 2, 1), more = FALSE)
```

<img src="figure/scaling-length-1.png" title="plot of chunk scaling-length" alt="plot of chunk scaling-length" width="700px" />

```r
#trellis.par.set("superpose.line", tpar$superpose.line)
```

## Multithreading efficiency

