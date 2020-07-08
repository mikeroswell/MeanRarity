# draft code to make plots that show mean rarity on a balanced "rarity scale." Goal is to make a flexible function that will take a vector of abundances and return a ggplot object that can be easily incorporated into figures.

# library(tidyverse)
# library(ggthemes)
# library(scales) # trans_new() is in the scales library

## TODO

## OPTIONAL: color the species (so that stacking becomes clear)

## get text and axis widths to scale with device/viewport. currently base_size argument can control manually

## think about squeezing when rarities aren't equal but are close s.t. boxes overlap... smart option that can calculate whether overlap occurs? Currently, have the lines=T option to plot species as lines if overlap is a problem (set by user)

#quick function to add together identical values

#' Sums duplicate values, for determining [species] weights
#'
#' Intended for use in constructing rarity plots, this function takes a numeric
#' vector of integers representing species abundances, and returns the weights
#' for each, unique rarity value. This is to solve overplotting and convey the
#' total "weight" of all species of a given rarity.
#'
#' @param x numeric vector
#'
#' @return A vector of weights, where only unique quantities have weights. The
#'   length of the new vector could be much shorter than \code{length(x)}.
#' @export
#' @examples
#' combfun(c(1,2,3,3,3,4,5, 9))
combfun<-function(x){x=x[order(x)]; y=unique(x)*as.numeric(table(x)); return (y)}


#' Round numbers, more aggressively the larger they are
#'
#' @param breaks is a numeric vector
#' intended for the breaks returned by \code{scales::transbreaks}
#'
#' @return numeric vector of same length as input, rounded nicely
#'
#' @export
#' @examples prettify(sqrt(2:7))
#' @examples prettify(sqrt(seq(2e3, 2e5, 4739)))
prettify <- function(breaks){
    digits <- -floor(log10(abs(breaks))) + 1
    digits[breaks == 0] <- 0
    raw<-round(breaks, digits = digits)
    raw[raw==-Inf]<-1e9
    # raw[is.na(raw)]<-1e11
    # raw[raw==0]<-1e-11
    return(raw) #klugey fix to -Inf
}

#' Power scale transformation.
#'
#' This is a function for nice scale transformations for \code{ggplot2},
#' including when \code{pow < 0}
#'
#' This function is for power transformations (i.e. raising to a power and
#' inverse) This is useful for visualizing generalized means.
#'
#' @param pow Exponent of power tranformation, scalar.
#' @param nb Number of desired breaks (approximate), scalar.
#'
#' @return A scale transformation object for plotting in ggplot.
#' @seealso \code{\link[scales]{trans_new}}, \code{\link[scales]{trans_breaks}},
#'   \code{\link{pfun}}, \code{\link{ipfun}}, \code{\link{prettify}}
#' @export
power_trans = function(pow, nb) scales::trans_new(name="power"
   , transform = function(x) pfun(x, pow)
   , inverse = function(x) ipfun(x, pow)
   , breaks = function(x) prettify(
		scales::trans_breaks(
			function(x) pfun(x, pow), function(x) ipfun(x, pow), n = nb
		)(x*1.1)
	)
   , domain = c(1, 10000) #this is to deal with -Inf
)

#' Replicate abundance vector to make stacks with geom_point
#'
#' @param df data.frame of abundances and rarities.
#'
#' @return data.frame with 1 row per individual for each species for plotting.
#'
#' @noRd
fancy_rep<-function(df){
    return(data.frame(df[rep(1:nrow(df), df$abundance),])
           %>% group_by(abundance)
           %>% mutate(gr=1:length(abundance), inds=rep(length(abundance)
               , length(abundance)))
    )
}

#' Base plot onto which balance plot is printed
#'
#' Takes abundance and various scale/dimension arguments to craft a ggplot
#' object for rarity plots
#'
#' @param abundance Numeric vector of (integer) species abundances.
#' @param fill_col Color for filling each of the boxes representing a single
#'   individual.
#' @param y_extent Scalar, how tall to draw y-axis.
#' @param x_max Scalar, approximately how far right should x-axis extend.
#' @param x_min Scalar, approximately how far left should x-axis extend.
#' @param base_size Typeface size for ggplot text (scalar).
#' @param noco Scalar, shrinks text and points if plotting muliple balance plots
#'   in a single plotting window
#' @param lines Logical, should each individual be plotted as a "box" or should
#'   individuals be summarized simply as the height of a line segment
#' @param verbose Logical, should the function return a pile of text
#'
#' @return ggplot object with some elements of a balance plot
#'
#' @noRd
base_plot <- function(abundance, pointScale
                      , fill_col="lightgrey" #can set to match communities
                      , y_extent=max(max(abundance),15) #how tall to draw y
                      , x_max=sum(abundance)/min(abundance) #plots a point to extend x
                      , x_min=sum(abundance)/max(abundance) # point to exted x
                      , base_size=24 #controls text size, default for 7" sq plotting device
                      , noco=1 #number of columns, shrinks text and point size proportional to number of panels
                      , lines=F
                      , verbose=T
                      ){
    #0.0353 is approximate points to cm conversion (a little less than 3 pts per mm)




    #make plotting data
    rf <- tibble(names = as.factor(1:length(abundance))
		, abundance
		, rarity = sum(abundance)/abundance
	)
	rfrepeated <-fancy_rep(rf)
	y_extent<-max(y_extent, max(combfun(abundance)))
	#14 is empirically derived scaling factor; but isn't quite right.
	# Seems like stuff below axis is about 2.5* height of 1 line of text
	pointScale<-(14*(min(dev.size("cm"))/noco-(2.5*0.0353*base_size)))
	pointsize <- pointScale/(y_extent*1.1)

	#0.5; shape is centered on  x,y; offset so it rests upon x, y-1
	goff <- 0.5

	#ggplot command to generate basic plot object
	base <- (ggplot(rfrepeated, aes(x=rarity, y=abundance))
	         +(if(lines==T){
	             rfdull<-rf %>% group_by(rarity) %>% summarize(inds=sum(abundance))
	   #line segments
	             geom_segment(data=rfdull, aes(x=rarity, xend=rarity, y=inds, yend=0)
	                           , color=rgb(0,0,0,0.4)
	                           , size=1
	                           )
	         } else{

	   #bricks
	    geom_point(aes(y=gr-goff, alpha=0.2), size=pointsize, fill=fill_col
		             , shape=22, color="black", stroke=0.5/noco)
	         })
		# plank
		+ geom_segment(
			aes(x, y, xend=xend, yend=yend)
			, data=data.frame(
				x=c(min(rf$rarity))
				, y=c(0)
				, xend=c(max(rf$rarity))
				, yend=c(0)
			)
		)

#fix plank location and add space above and below data range
		+scale_y_continuous(
		    expand=c(0,0)
		    , limits=c(y_extent-1.05*y_extent, 1.05*y_extent)
		)
    	+ labs(y="individuals")
	)
	#calls the function theme plot to generate basic figure
	return(theme_plot(base, base_size=base_size, noco=noco))
}

#' Set preferences for axes etc. for balance plots
#'
#' @param p ggplot object
#' @param base_size Typeface size for ggplot text (scalar).
#' @param noco Scalar, shrinks text and points if plotting muliple balance plots
#'   in a single plotting window
#'
#' @return also a ggplot object, with some theme elements specified
#'
#' @noRd
theme_plot <- function(p, base_size=24, noco=1,...){
	return(p
		+ theme_tufte(base_family = "sans", base_size=base_size/noco)
		+ theme(legend.position="none")
		+ theme(
		    axis.text=element_text(color="black") # better than Wickham's grey?
		    , axis.line.x = element_line(colour = 'black'
		                                 , size=0.2, linetype='solid')
			, axis.line.y = element_line(colour = 'black'
			                             , size=0.2, linetype='solid')
		)
	)
}

#' Rescales baseplot x-axis, adds space to allow matching scales between comms
#'
#' @param ell exponent for scaling generalized mean
#' @param fill_col Color for filling each of the boxes representing a single
#'   individual.
#' @param y_extent Scalar, how tall to draw y-axis.
#' @param x_max Scalar, approximately how far right should x-axis extend.
#' @param x_min Scalar, approximately how far left should x-axis extend.
#' @param noco Scalar, shrinks text and points if plotting muliple balance plots
#'   in a single plotting window
#' @param lines Logical, should each individual be plotted as a "box" or should
#'   individuals be summarized simply as the height of a line segment
#' @param nbreaks Integer, approximate number of x-axis tick marks
#'
#' @return a gpplot object with some of the theme items set
#'
#' @noRd
scale_plot <- function(
        ab
        , ell
        , fill_col="lightgrey"
        , y_extent=max(max(ab), 15)
        , x_max=sum(ab)/min(ab)
        , x_min=sum(ab)/max(ab)
        , noco=1
        , lines=F
        , nbreaks=5, ...
){
    return (base_plot(ab, fill_col=fill_col, y_extent=y_extent
	            , x_max=x_max, x_min=x_min, noco=noco, lines=lines, ...)
		+ scale_x_continuous(trans=power_trans(pow=ell, nb=nbreaks), labels=signif)
		+ geom_point(aes(x,y) #allows for x min and max points to determine axes
		  , data=tibble(x=c(x_max, x_min), y=c(0,0))
		  , color="white", alpha=0)
	        )
}

#' Reference points at means with power \code{ell}
#'
#' @param ab Numeric vector of integer species abundances
#' @param ell Exponent for type of mean rarity (scalar)
#' @param noco Scalar, shrinks text and points if plotting muliple balance plots
#'   in a single plotting window
#'
#' @return geom object to add to a ggplot object to construct balance plot
#'
#' @noRd
mean_points <- function(ab, ell, noco=1){
    ab<-ab[ab!=0]
	div <- Vectorize(dfun, vectorize.args=("l"))(ab, ell)
	return(geom_point(
		data=tibble(x=div, y=0*div, clr=1:length(div))
		, aes(x, y, color=as.factor(clr))
		, size=0.2*min(dev.size("cm"))/noco
	))
}

#' Plot the fulcrum
#'
#' @noRd
fulcrum<-function(ab, ell, y_extent=max(max(combfun(ab)), 15), x_max=1
                  , x_min=1, fill_col="light_grey"
                  , base_size=24, noco=1, nbreaks=5, verbose=T){

    ab<-ab[ab!=0]
    div <- dfun(ab, ell)

    if(verbose==T){
        print(c(paste("diversity =", div), paste("community size =", sum(ab))
                , paste("max observed rarity =", sum(ab)/min(ab))
                , paste("min observed rarity =", sum(ab)/max(ab))))}

    return(geom_point(
        data=tibble(x=div, y=-0.03*y_extent) # gets fulcrum point close.
        , size=(0.48*min(dev.size("cm"))-(2.5*0.0353*base_size))/noco #scales with plotting device and number of columns
        # , size=rel(0.3)
        , shape=17
        , aes(x, y)
    )
    )
}

#' Construct the full balance plot for scale ell, with reference means=means
#'
#' @seealso This function depends on internal functions in the \code{MeanRarity}
#'   package which can be accessed with \code{:::} e.g.
#'   \code{MeanRarity:::scale_plot}
#'
#'
#' @export


rarity_plot <- function(ab
                        , ell
                        , means=-1:1
                        , noco=1
                        , lines=F
                        , ...){
    ab<-ab[ab!=0]
    print(cat("     rarity plot expects a square viewport and resizes points based on\n     min(dev.size() and noco (for number of columns).\n     Selecting lines=T will plot stacks of individuals as a line element,\n     which tends to be more robust to window size.\n     lines=T may be the best way to deal with overplotting,\n     which results from several species with similar but not identical rarities.\n "))
	return(
		scale_plot(ab
		           , ell
		           , noco=noco
		           , lines=lines
		           ,...)
		+ mean_points(ab
		              , means
		              , noco=noco)
		+ fulcrum(ab
		          , ell
		          , noco=noco
		          , ...)
		+ scale_color_brewer(type="qual", palette="Set1")
		# + scale_color_viridis_d(option="plasma")
		# +scale_color_manual(values=c("#AA8E39", "#294F6D", "#4B2D73"))
	)
}

#conventiently plot for l=-1:1, with reference points in each fig
rarity_series <- function(ab, lrange=-1:1, means=lrange,...){
	for(l in lrange){
		print(rarity_plot(ab, l, means,...))
	}
}

#convenience functions to omit y-axis elements for constructing multi-panel plots.
white_y<-function(p){
    return(p
           +theme(axis.text.y=element_text(color="white")
                  , axis.title.y=element_text(color="white")
                  , axis.ticks.y = element_line(color="white")
                  , axis.line.y = element_line(color="white")
                  , axis.line.x=element_line(
                      colour = 'black', size=0.2, linetype='solid'
                    )
                  )
           )
}

omit_y<-function(p){
    return(p
           +theme(axis.text.y = element_blank()
                  , axis.title.y = element_blank()
                  , axis.ticks.y = element_blank()
                  , axis.line.y = element_blank()
                  , axis.line.x=element_line(
                      colour = 'black', size=0.2, linetype='solid'
                  )
           )
    )
}

#some SADs to play with


# ab <- c(20, 15, 9, 3, 2, 1, 1) #includes stacking
ab<-c(20,8,5,4,2,1) #candidate for user's guide
# ab <- c(100, 20, 15, 9, 3, 2, 1, 1)
# ab <- c(50,30,20,0,0,0)
# ab <- c(4,3,2)
# ab <- c(20, 15, 9, 3, 2, 1, 1, 0, 0)
# ab <- c(200,100, 20, 15, 9, 3, 2, 1, 1)
# ab <- floor(exp(rnorm(50, 4, 1.5)))


#################
# code for figs in users guide etc.


# rarity_plot(ab,0, fill_col="red", base_size=24, verbose=T, noco=2)
#
# p<-rarity_plot(ab, 1, fill_col="blue", x_min=1, x_max=45, noco=3, base_size=12)
#
# grid.arrange(p, omit_y(p), omit_y(p), p, omit_y(p), omit_y(p), p, omit_y(p), omit_y(p))

# rarity_series(ab=ab, 1:-1)
#
# plot some rarity "balance plots" for guide to measuring diversity (Roswell et al. 2020 Oikos....)
# pdf(file="figures/rarity_plots_for_guide.pdf")
# rarity_plot(ab, 1)+scale_color_brewer(type="qual", palette="Dark2") #also changed base_size to 30 or maybe 32 and doubled point size.
# white_y(rarity_plot(ab,0))+scale_color_brewer(type="qual", palette="Dark2")
# white_y(rarity_plot(ab,-1))+scale_color_brewer(type="qual", palette="Dark2")
# dev.off()
##############
# RAD plot code
radplot<-function(comm
                  , maxrich=length(comm)
                  , maxab=max(comm)
                  , fill
                  , shape=16){
    comm<-comm[comm!=0]
    rawrnk<-tibble(abund=comm, rnk=row_number(comm))
    toplot<-rawrnk %>%
        mutate(x=-rnk-maxrich+max(rnk))

    f<-(toplot %>% ggplot(aes(x, abund, size))
        +geom_point(shape=shape, color=fill, size=2)
        +geom_line(color=fill)
        + scale_x_continuous(limits=c(-maxrich, 0))
        + scale_y_continuous(limits=c(0,maxab))
        +theme_classic()
        +ggtitle("              rank-abundance plot")
        +theme(axis.text.x = element_text(color="white"), axis.text.y=element_text(colour="black")
               , legend.position="none", text=element_text(size=12))
        +labs(x="abundance rank", y="individuals")
    )
    return(f)

}
