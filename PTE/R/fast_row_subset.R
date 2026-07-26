#' Subset rows (and optionally columns) of a data frame without \code{[.data.frame}'s overhead
#'
#' \code{Xy[idx, ]} (and especially \code{Xy[sample(1:n, n, replace = TRUE), ]}) is expensive
#' when called tens of thousands of times across the bootstrap x cross-validation loop: it
#' re-derives factor levels, re-dispatches through \code{[[.data.frame}/\code{[.data.frame} for
#' every column, and -- worst of all when there are repeated indices from bootstrap resampling --
#' calls \code{make.unique()} to disambiguate duplicated row names via \code{paste}/\code{unique}.
#' None of that bookkeeping is needed here since nothing downstream looks up results by row name.
#' This does the same per-column extraction (still dispatching to \code{[.factor} etc. so factor
#' levels are preserved correctly) but skips the data-frame-level machinery, and is safe to use
#' regardless of which \code{personalized_model_build_function}/\code{predict_function} the caller
#' plugs in since the result is a normal, valid data frame (just with plain integer row names).
#'
#' @param df 	A data frame.
#' @param idx 	A row index vector (positive or negative, exactly as you'd pass to \code{df[idx, ]}).
#' @param cols 	A column index vector (defaults to all columns).
fast_row_subset = function(df, idx, cols = seq_along(df)){
	n_out = length(attr(df, "row.names")[idx])
	out = vector("list", length(cols))
	for (k in seq_along(cols)){
		col = .subset2(df, cols[k])
		#a data frame column can itself be an embedded matrix (e.g. from scale() or poly());
		#plain vector-style col[idx] would silently flatten/corrupt that, so subset rows only
		out[[k]] = if (is.matrix(col)) col[idx, , drop = FALSE] else col[idx]
	}
	names(out) = names(df)[cols]
	attr(out, "row.names") = seq_len(n_out)
	class(out) = "data.frame"
	out
}
