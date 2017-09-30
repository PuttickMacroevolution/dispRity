    
#     if(plot.disparity) {
    
#         par(mfrow=c(1, 2), mar=c(4,4,2,2), oma=c(10, 4, 4, 4))
#         xAxis <- max(model_test_input[[4]]) - model_test_input[[4]]
#         plot(xAxis, model_test_input[[1]], type="l", xlim=c(max(xAxis), 0), xlab="Time", ylab="central tendency", las=1)
#         varUp <- model_test_input[[1]] + model_test_input[[2]]
#         varDown <- model_test_input[[1]] - model_test_input[[2]]
#         polygon(x=c(xAxis, rev(xAxis)), c(varUp, rev(varDown)), col="grey", border=F)
#         lines(xAxis, model_test_input[[1]], col="grey50")
#         abline(v=time.split, lty=2, lwd=2)
#         plotcor <- barplot(weight.aicc[order.aicc], las=1, ylim=c(0, 1), col="grey30", border=F, ylab="Akaike weights", names=F)
#         mtext(model.names[order.aicc], 1, las=2, at=plotcor[,1], line=1)
#     }
    

    
# plot.disparity.time <- function(input_disparity, plot.variance=FALSE, plot.coords=FALSE) {
	
# 	model.test_input <- select.model.list(input_disparity)
# 	xAxis <- max(model.test_input[[4]]) - model.test_input[[4]]
# 	plot(xAxis, model.test_input[[1]], type="l", xlim=c(max(xAxis), 0), xlab="Time", ylab="disparity measure", las=1)
# 	if(plot.variance) {
# 		varUp <- model.test_input[[1]] + model.test_input[[2]]
# 		varDown <- model.test_input[[1]] - model.test_input[[2]]
# 		polygon(x=c(xAxis, rev(xAxis)), c(varUp, rev(varDown)), col="grey", border=F)
# 	}	
# 	lines(xAxis, model.test_input[[1]], col="grey50")
# 	if(plot.coords) return(cbind(xAxis, model.test_input[[1]]))
# }

	
# 	