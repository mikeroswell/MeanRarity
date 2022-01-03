# draft code to make plots that show mean rarity on a balanced "rarity scale."
# Goal is to make a flexible function that will take a vector of abundances and
# return a ggplot object that can be easily incorporated into figures.

# library(tidyverse)
# library(ggthemes)
# library(scales) # trans_new() is in the scales library

## TODO

## OPTIONAL: color the species (so that stacking becomes clear)

## get text and axis widths to scale with device/viewport. currently base_size
## argument can control manually

## think about squeezing when rarities aren't equal but are close s.t. boxes
## overlap... smart option that can calculate whether overlap occurs? Currently,
## have the lines=T option to plot species as lines if overlap is a problem (set
## by user)


#' Sums duplicate values, for determining \[species\] weights
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
#' combfun(c(1, 2, 3, 3, 3, 4, 5, 9))
combfun <- function(x){
  x = x[order(x)]
  y = unique(x) * as.numeric(table(x))
  return(y)}


#' Round numbers, more aggressively the larger they are
#'
#' @param breaks is a numeric vector
#' intended for the breaks returned by \code{scales::transbreaks}
#'
#' @return numeric vector of same length as input, rounded nicely
#'
#' @noRd
#' @examples prettify(sqrt(2:7))
#' @examples prettify(sqrt(seq(2e3, 2e5, 4739)))
#'
prettify <- function(breaks){
    digits <- -floor(log10(abs(breaks))) + 1
    digits[breaks == 0] <- 0
    raw <- round(breaks, digits = digits)
    raw[raw == -Inf] <- 1e9
    # raw[is.na(raw)]<-1e11
    # raw[raw==0]<-1e-11
    return(raw) #klugey fix to -Inf
}

#' Power scale transformation
#'
#' This is a function for nice scale transformations for \code{ggplot2},
#' including when \code{pow < 0}.
#'
#' This function is for power transformations (i.e. raising to a power and
#' inverse) This is useful for visualizing generalized means.
#'
#' @param pow Exponent of power transformation, scalar.
#' @param nb Number of desired breaks (approximate), scalar.
#'
#' @return A scale transformation object for plotting in ggplot.
#' @seealso \code{\link[scales]{trans_new}}, \code{\link[scales]{trans_breaks}},
#'   \code{\link{pfun}}, \code{\link{ipfun}}
#' @export
power_trans = function(pow, nb) scales::trans_new(name = "power"
   , transform = function(x) pfun(x, pow)
   , inverse = function(x) ipfun(x, pow)
   , breaks = function(x) prettify(
		scales::trans_breaks(
			function(x) pfun(x, pow), function(x) ipfun(x, pow), n = nb
		)(x * 1.1)
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
fancy_rep <- function(df){
  mydf = data.frame(df[rep(1:nrow(df), df$ab), ])
    return(mydf %>%
              dplyr::group_by(.data$ab) %>%
              dplyr::mutate(gr = 1:length(.data$ab)
                            , inds = rep(length(.data$ab), length(.data$ab)))
    )
}

#' Base plot onto which balance plot is printed
#'
#' Takes abundance and various scale/dimension arguments to craft a ggplot
#' object for rarity plots.
#'
#' @template ab_template
#' @param fill_col Color for filling each of the boxes representing a single
#'   individual.
#' @param y_extent Scalar, how tall to draw y-axis.
#' @param x_max Scalar, approximately how far right should x-axis extend.
#' @param x_min Scalar, approximately how far left should x-axis extend.
#' @param base_size Typeface size for ggplot text (scalar).
#' @param noco Scalar, shrinks text and points if plotting multiple balance
#'   plots in a single plotting window
#' @param lines Logical, should each individual be plotted as a "box" or should
#'   individuals be summarized simply as the height of a line segment
#' @param verbose Logical, should the function return a pile of text
#'
#' @return ggplot object with some elements of a balance plot
#'
#' @noRd
base_plot <- function(ab, pointScale
                      , fill_col = "lightgrey" #can set to match communities
                      , y_extent = max(max(ab), 15) #how tall to draw y
                      , x_max = sum(ab)/min(ab) #plots a point to extend x
                      , x_min = sum(ab)/max(ab) # point to extend x
                      , base_size = 24 #controls text size, default for 7" sq plotting device
                      , noco = 1 #number of columns, shrinks text and point size proportional to number of panels
                      , lines = F
                      , verbose = T
                      ){

    #make plotting data
    rf <- tibble::tibble(names = as.factor(1:length(ab))
		, ab
		, rarities = sum(ab)/ab
	)
	rfrepeated <- fancy_rep(rf)
	y_extent <- max(y_extent, max(combfun(ab)))
	#14 is empirically derived scaling factor; but isn't quite right.
	# Seems like stuff below axis is about 2.5* height of 1 line of text
	typeface_size = 14

	pointScale <- (typeface_size *
	                 (min(grDevices::dev.size("cm")) /
	                    noco - (2.5 * 1/ggplot2::.pt/10 * base_size)
	                  )
	               )
	pointsize <- pointScale / (y_extent * 1.1)

	#0.5; shape is centered on  x,y; offset so it rests upon x, y-1
	goff <- 0.5

	#ggplot command to generate basic plot object
	base <- (rfrepeated %>%
	           ggplot2::ggplot(ggplot2::aes(x = .data$rarities, y = ab)) +
	         (if(lines == T){
	             rfdull <- rf %>%
	               dplyr::group_by(.data$rarities) %>%
	               dplyr::summarize(inds = sum(.data$ab))
	   #line segments
	             ggplot2::geom_segment(data = rfdull,
	                                   ggplot2::aes(x = .data$rarities
                                                  , xend = .data$rarities
                                                  , y = .data$inds
	                                                , yend = 0)
	                           , color = grDevices::rgb(0,0,0,0.4)
	                           , size = 1
	                           )
	         } else{

	   #bricks
	           ggplot2::geom_point(ggplot2::aes(y = .data$gr - goff, alpha = 0.2)
	                               , size = pointsize
	                               , fill = fill_col
		                             , shape = 22
		                             , color = "black"
		                             , stroke = 0.5/noco)
	         })
		# plank
		+ ggplot2::geom_segment(
		  ggplot2::aes(x = .data$x, y = .data$y
		               , xend = .data$xend, yend = .data$yend)
			, data = data.frame(
				x = c(min(rf$rarities))
				, y = c(0)
				, xend = c(max(rf$rarities))
				, yend = c(0)
			)
		)

#fix plank location and add space above and below data range
		+ggplot2::scale_y_continuous(
		    expand = c(0,0)
		    , limits = c(y_extent-1.05 * y_extent, 1.05 * y_extent)
		)
    	+ ggplot2::labs(y = "individuals")
	)
	#calls the function theme plot to generate basic figure
	return(theme_plot(base, base_size = base_size, noco = noco))
}

#' Set preferences for axes etc. for balance plots
#'
#' @param p ggplot object
#' @param base_size Typeface size for ggplot text (scalar).
#' @param noco Scalar, shrinks text and points if plotting multiple balance plots
#'   in a single plotting window
#' @param ... Additional arguments passed to other functions.
#'
#' @return A ggplot object, with some theme elements specified
#'
#' @noRd
theme_plot <- function(p, base_size = 24, noco = 1, ...){
	return(p
		+ ggthemes::theme_tufte(base_family = "sans", base_size=base_size/noco)
		+ ggplot2::theme(legend.position = "none")
		+ ggplot2::theme(
		    axis.text = ggplot2::element_text(color = "black") # better than Wickham's grey?
		    , axis.line.x = ggplot2::element_line(colour = 'black'
		                                 , size = 0.2, linetype = 'solid')
			, axis.line.y = ggplot2::element_line(colour = 'black'
			                             , size = 0.2, linetype = 'solid')
		)
	)
}

#' Rescales baseplot x-axis, adds space to allow matching scales between comms
#'
#' @template l_template
#' @param fill_col Color for filling each of the boxes representing a single
#'   individual.
#' @param y_extent Scalar, how tall to draw y-axis.
#' @param x_max Scalar, approximately how far right should x-axis extend.
#' @param x_min Scalar, approximately how far left should x-axis extend.
#' @param noco Scalar, shrinks text and points if plotting multiple balance
#'   plots in a single plotting window
#' @param lines Logical, should each individual be plotted as a "box" or should
#'   individuals be summarized simply as the height of a line segment
#' @param nbreaks Integer, approximate number of x-axis tick marks
#' @param ... Additional arguments passed to other functions.
#'
#' @return a gpplot object with some of the theme items set
#'
#' @noRd
scale_plot <- function(
        ab
        , l
        , fill_col = "lightgrey"
        , y_extent = max(max(ab), 15)
        , x_max = sum(ab)/min(ab)
        , x_min = sum(ab)/max(ab)
        , noco = 1
        , lines = F
        , nbreaks = 5
        , ...
){
    return(base_plot(ab
                      , fill_col = fill_col
                      , y_extent = y_extent
	                    , x_max = x_max
                      , x_min = x_min
                      , noco = noco
                      , lines = lines
                      , ...)
		+ ggplot2::scale_x_continuous(trans = power_trans(
		    pow = l, nb = nbreaks), labels = signif)
		+ ggplot2::geom_point(ggplot2::aes(.data$x, .data$y) #allows for x min and max points to determine axes
		  , data = tibble::tibble(x = c(x_max, x_min), y = c(0,0))
		  , color = "white", alpha=0)
	        )
}

#' Reference points at means with power \code{l}
#'
#' @template ab_template
#' @template l_template
#' @param noco Scalar, shrinks text and points if plotting multiple balance
#'   plots in a single plotting window.
#'
#' @return geom object to add to a ggplot object to construct balance plot
#'
#' @noRd
mean_points <- function(ab, l, noco = 1){
    ab <- ab[ab != 0]
	div <- Vectorize(rarity, vectorize.args = ("l"))(ab, l)
	pointDat = tibble::tibble(x = div, y = 0*div, clr = 1:length(div))
	return(ggplot2::geom_point(data = pointDat
	  , ggplot2::aes(x = .data$x, y = .data$y, color = as.factor(.data$clr))
		, size = 0.2 * min(grDevices::dev.size("cm")) / noco
	))
}

#' Plot the fulcrum
#'
#' @template l_template
#' @param fill_col Color for filling each of the boxes representing a single
#'   individual.
#' @param y_extent Scalar, how tall to draw y-axis.
#' @param x_max Scalar, approximately how far right should x-axis extend.
#' @param x_min Scalar, approximately how far left should x-axis extend.
#' @param noco Scalar, shrinks text and points if plotting multiple balance
#'   plots in a single plotting window.
#' @param nbreaks Integer, approximate number of x-axis tick marks.
#' @param verbose Logical, should the function return a pile of text.
#'
#'
#' @noRd
fulcrum <- function(ab, l
                  , y_extent = max(max(combfun(ab)), 15)
                  , x_max = 1
                  , x_min = 1
                  , fill_col = "light_grey"
                  , base_size = 24
                  , noco = 1
                  , nbreaks = 5
                  , verbose = T){

    ab <- ab[ab != 0]
    div <- rarity(ab, l)

    if(verbose == T){
        print(c(paste("diversity =", div), paste("community size =", sum(ab))
                , paste("max observed rarity =", sum(ab) / min(ab))
                , paste("min observed rarity =", sum(ab) / max(ab))))}
    # y_off = -0.03 # small amount to offset the y axis.
    y_off = -1/ggplot2::.pt/10
    size_adjust = 0.5 # this shrinks the fulcrum given the size of other
    # points for shape = 17 to have its apex at y = 0.
    fulcDat = tibble::tibble(x = div, y = y_off * y_extent)
    return( ggplot2::geom_point(data = fulcDat
        , size = (size_adjust * min(grDevices::dev.size("cm")) -
                    (2.5 * 1/ggplot2::.pt/10 * base_size)) /
            noco # scales with plotting device and number of columns
        # , size=rel(0.3)
        , shape = 17
        , ggplot2::aes(.data$x, .data$y)
    )
    )
}

#' Construct rarity balance plot
#'
#' This function takes the abundance vector, scaling exponent, and target means
#' (Default the Pythagorean means), and returns a formatted 1-panel ggplot
#' object.
#'
#' Hill diversity, or "mean rarity," is the balance point for the community
#' along the rarity scale. The image produced by \code{rarity_plot} illustrates
#' this balance. Each block represents an individual: because Hill
#' diversities are weighted by abundance, the “mass” of each “block” is the same
#' regardless of species identity. Each individual’s x-axis value is given its
#' species’s "rarity," which is the reciprocal of its relative abundance. The
#' parameter \code{l} controls how rarity is scaled. A community’s balance
#' point along the rarity scale, pictured as a triangular fulcrum, is the mean
#' rarity, or diversity, of the community.
#'
#' To ease comparison across scales, by default the Pythagorean means
#' (\url{https://en.wikipedia.org/wiki/Mean}) are marked
#' with reference points: the arithmetic mean with a rose dot, the geometric
#' mean with a blue dot, and the harmonic mean with a green dot. The arithmetic
#' scale provides high leverage to very rare species; although they carry little
#' weight (few individuals), these species influence the mean a great deal
#' because they sit far to the right of the rarity scale. The arithmetic mean
#' rarity of the community is the Hill diversity when \eqn{\ell = 1}{ℓ = 1}, and is equal to
#' species richness. The logarithmic scale provides less leverage to very rare
#' species. Thus, the geometric mean rarity of the community is lower. The
#' geometric mean rarity is also known as the Hill-Shannon diversity, or the
#' Hill diversity when \eqn{\ell = 0}{ℓ = 0}. The reciprocal scale accords more leverage to low
#' rarity values. Thus, the harmonic mean rarity, also known as the Hill-Simpson
#' diversity, or Hill diversity when \eqn{\ell = -1}{ℓ = -1}, is much lower still. An
#' interactive online application that enables users to specify species
#' abundances and the scaling parameter is available at
#' \url{https://mean-rarity.shinyapps.io/rshiny_app1/}
#'
#' The scaling of various plot elements depends on the plotting device.
#' Rstudio's seems especially touchy based on window size. The defaults here
#' pertain to the standard 7"x7" plotting window given by \code{quartz()} or
#' \code{pdf()}. They also seem to play nice on the shiny app
#' \url{https://mean-rarity.shinyapps.io/rshiny_app1/}. Other window sizes or
#' devices may require tweaking.
#'
#' @template ab_template
#' @template l_template
#' @template q_template
#' @param noco Scalar, shrinks text and points if plotting multiple balance
#'   plots in a single plotting window.
#' @param lines Logical, should each individual be plotted as a "box" or should
#'   individuals be summarized simply as the height of a line segment.
#' @param means Numeric vector of scaling exponent values corresponding to
#'   reference (by default, Pythagorean) means.
#' @param ... Additional arguments passed to other functions.
#'
#'
#' @seealso This function depends on internal functions in the
#'   \code{\link{MeanRarity}} package which can be accessed with \code{:::} e.g.
#'   \code{MeanRarity:::scale_plot}.
#'
#' @concept Visualization
#'
#' @export
#' @examples
#' ab<-c(20,8,5,4,2,1)
#'
#' # experiment with other abundance vectors!
#' # ab <- c(20, 15, 9, 3, 2, 1, 1)
#' # ab <- c(100, 20, 15, 9, 3, 2, 1, 1)
#' # ab <- c(50,30,20,0,0,0)
#' # ab <- c(4,3,2)
#' # ab <- c(20, 15, 9, 3, 2, 1, 1, 0, 0)
#' # ab <- c(200,100, 20, 15, 9, 3, 2, 1, 1)
#' # ab <- floor(exp(rnorm(50, 4, 1.5)))
#'
#' richness <- rarity_plot(ab, 1)
#' Hill_Shannon <- rarity_plot(ab, 0)
#' Hill_Simpson <- rarity_plot(ab, -1)
#'
#' richness
#' Hill_Shannon
#' Hill_Simpson
#'
#' # richness + Hill_Shannon + Hill_Simpson # plot with patchwork
rarity_plot <- function(ab
                        , l
                        , q = NULL
                        , means = -1:1
                        , noco = 1
                        , lines = FALSE
                        , ...){
 message(cat("     rarity plot expects a square viewport (likely issues in the RStudio plotting device) and resizes points based on\n     min(dev.size() and noco (for number of columns).\n     Selecting lines = T will plot stacks of individuals as a line element,\n     which tends to be more robust to window size.\n     Setting lines = T may be the best way to deal with overplotting,\n     which results from several species with similar but not identical rarities. \n "))
  if(!is.null(q)){
    l = 1-q
    warning("l has been set to 1-q")
  }

    ab <- ab[ab != 0]
	return(
		scale_plot(ab
		           , l
		           , noco=  noco
		           , lines = lines
		           ,...)
		+ mean_points(ab
		              , means
		              , noco = noco)
		+ fulcrum(ab
		          , l
		          , noco = noco
		          , ...)
		+ ggplot2::scale_color_brewer(type = "qual", palette = "Set1")
		# + scale_color_viridis_d(option="plasma")
		# + scale_color_manual(values=c("#AA8E39", "#294F6D", "#4B2D73"))
	)
}

#' Conveniently plot for l = -1:1, with reference points in each fig
#'
#' Convenience function for plotting balance plots for the arithmetic,
#' geometric, and harmonic mean rarities with only one line of code.
#'
#' @return Prints 3 ggplots into graphics device.
#'
#' @template l_template
#' @param lrange Numeric vector of scaling exponent values
#' @param ... Additional arguments passed to other functions.
#'
#' @noRd
#'
rarity_series <- function(ab, lrange = -1:1, means = lrange, ...){
	for(l in lrange){
		print(rarity_plot(ab, l, means, ...))
	}
}

#' Convenience function to blank y-axis elements for constructing multi-panel
#' plots
#'
#' This function whites out, but does not remove, y-axis elements for a
#' consistent plot size that can be included in a multi-panel plot
#'
#' @return A ggplot.
#'
#' @param p A ggplot.
#'
#' @export
white_y<-function(p){
    return(p
           + ggplot2::theme(axis.text.y = ggplot2::element_text(color="white")
                  , axis.title.y = ggplot2::element_text(color="white")
                  , axis.ticks.y = ggplot2::element_line(color="white")
                  , axis.line.y = ggplot2::element_line(color="white")
                  , axis.line.x = ggplot2::element_line(
                      colour = 'black', size=0.5, linetype='solid'
                    )
                  )
           )
}


#' Convenience function to omit y-axis elements for constructing multi-panel
#' plots
#'
#' This function removes y-axis elements for plots intended to be
#' included in a multi-panel plot. When the y-axis is omitted (rather than
#' whited-out as in \code{white_y}), the plotted area typically expands to fill
#' unused space.
#'
#' @seealso \code{\link{white_y}}
#'
#' @return A ggplot.
#'
#' @param p A ggplot.
#'
#' @export
omit_y<-function(p){
    return(p
           + ggplot2::theme(axis.text.y = ggplot2::element_blank()
                  , axis.title.y = ggplot2::element_blank()
                  , axis.ticks.y = ggplot2::element_blank()
                  , axis.line.y = ggplot2::element_blank()
                  , axis.line.x = ggplot2::element_line(
                      colour = 'black', size=0.2, linetype='solid'
                  )
           )
    )
}






#' Plot a rank-abundance distribution
#'
#' Take a vector of species \[relative\] abundances and plot a rank-abundance
#' distribution or Whittaker plot.
#'
#' @template ab_template
#' @param maxrich Scalar, how many species to include space for
#' @param maxab Scalar, how big is the largest abundance value
#' @param fill Character string specifying a color, could be extended to
#' permit name for community ID variable
#' @param shape Character string for point shape, could be extended to
#' permit name for community ID variable
#' @param size Scalar, adjusts point size in ggplot object
#' @param Whittaker Logical, if \code{TRUE} log-transform abundances
#'
#' @return ggplot object with a rank-abundance plot
#'
#' @concept Visualization
#' @export
#' @examples radplot(c(20,8,5,4,2,1))
#' radplot(c(20,8,5,4,2,1), Whittaker = TRUE)

radplot <- function(ab
                  , maxrich = length(comm)
                  , maxab = max(comm)
                  , fill = "red"
                  , shape = 16
                  , size = 2.5
                  , Whittaker = FALSE
                  ){
    comm = ab[ab!=0]
    if(Whittaker){comm = log(comm)}
    rawrnk = tibble::tibble(abund = comm, rnk = dplyr::row_number(comm))
    toplot = rawrnk %>%
        dplyr::mutate(x = -.data$rnk - maxrich + max(.data$rnk))

    f <- (toplot %>% ggplot2::ggplot(ggplot2::aes(.data$x, .data$abund, .data$size))
        + ggplot2::geom_point(shape = shape, color = fill, size = size)
        + ggplot2::geom_line(color = fill)
        + ggplot2::scale_x_continuous(limits = c(-maxrich, 0))
        + ggplot2::scale_y_continuous(limits = c(0, maxab))
        + ggplot2::theme_classic()
        + ggplot2::ggtitle(ifelse(Whittaker
                                  , "              rank-log(abundance) [Whittaker] plot"
                                  , "              rank-abundance plot"
                                  )
        )
        + ggplot2::theme(axis.text.x = ggplot2::element_text(color = "white")
                         , axis.text.y = ggplot2::element_text(colour = "black")
                         , legend.position = "none"
                         , text = ggplot2::element_text(size = 12)
                         )
        + ggplot2::labs(x = "abundance rank", y = ifelse(Whittaker
                                                         , "ln(individuals)"
                                                         , "individuals")
        )
    )
    return(f)

}
