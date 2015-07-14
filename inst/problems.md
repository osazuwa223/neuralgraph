Original process:
1. Make a DAG with hidden variables
2. Fit weights many times
3. Get an average of betweenness rankings
4. Prioritize hidden nodes by there betweeness
5. Return the hidden node and its MB

Problems:
	I. Optimization is not working.  Possible approach:
		1.  Collapse out the hidden nodes
		2.  Fit weights using NLS.  
			1. Use observed values to get a set of initial weights
			2. Propagate predictions
			3. Iteratively fit on last iterations predictive values, where each iteration has an updated penalty term 
		Notes: I see no simple way of making thiw work with hidden nodes.
		
	II. Betweenness is less meaningful when weights are not identifiable.  
		Possible approach1
			1. Collapse out the hidden nodes.  You planned to do this anyway, now you are just combining several 
			hidden nodes until 2
		Possible approach2
			1. Choose based on average ranking.  This assumes that the average ranking is stabalizes, and this
			depends on the nature of the search space.  
		Possible approach3.  
			1. Use predictive ability instead, and return markov blankets of observed.  See problem 4.



	III. If I collapse: Before I was chosing hidden nodes that had the highest betweenness.  If I revert to collapsing the graph, how will I select? 
		Approach 1
			1. Select the markov neighborhood of any observed node.  Quantify those.
	IV. If I collapse: before I was doing betweeness on uncollapsed network, how do I do betweeness on the collapse?
		Approach 1
			1.  Choose betweenness based on the collapsed network
	V. What looks like a better package?
		1.  Biologists would want to work with the original network in tact.  Creating a collapsed "dwarf" network
			would be unappealing.



What is not a problem
	I. Using coefficients as weights makes sense, as they represent effect size, as well as amount of information, if not type of informaiton.
	II.  I was thinking that prediction depends on the number of hidden nodes.  This is true, but I am comparing two models for selection based on the same input graph, so it doesn't matter.  Plus, I am penalizing.