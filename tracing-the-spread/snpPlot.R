data = read.csv(file = "entropies.csv")
xx <- barplot(height = data$entropy[data$entropy > 0])
axis(1, at=xx, labels=data$location[data$entropy > 0], tick=FALSE, las=2, line=-0.5, cex.axis=1.5)
text(x = xx, y = data$entropy[data$entropy > 0], label = round(data$entropy[data$entropy > 0], 3), pos = 3, cex = 0.8, col = "blue")

xx <- barplot(height = data$entropy)
axis(1, at=xx, labels=data$location, tick=FALSE, las=2, line=-0.5, cex.axis=0.75)
text(x = xx, y = data$entropy, label = round(data$entropy, 2), pos = 3, cex = 0.3, col = "blue")
