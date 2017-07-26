
set SERVERS;
set PAIRS;
set LCELLS;
set RCELLS;

param lcell_network{SERVERS, LCELLS};
param rcell_network{SERVERS, RCELLS};
param pair_weight{PAIRS};
param pair_left{PAIRS};
param pair_right{PAIRS, RCELLS};
param priour_mkspn{SERVERS};
param f;
set prs dimen 2 := setof{p in PAIRS, r in RCELLS: pair_right[p,r] = 1}(p,r);

var map{PAIRS,SERVERS}, binary; # >= 0;
var use_r{RCELLS, SERVERS} >= 0;
var use_l{LCELLS, SERVERS} >= 0;
var x0;

minimize z:
	/* communication */
	  (sum{r in RCELLS, s in SERVERS} (use_r[r,s] * rcell_network[s,r]) )
	+ (sum{l in LCELLS, s in SERVERS} (use_l[l,s] * lcell_network[s,l]) )
	
	/* balance */
	+ f * x0 
	;

#s.t. min_x0: x0 <=  402558;
s.t. only_one_server{p in PAIRS}: sum{s in SERVERS} map[p,s] = 1;
s.t. find_max{s in SERVERS}: x0 >= priour_mkspn[s] + sum{p in PAIRS} (map[p,s] * pair_weight[p]);
s.t. use_l_cells{s in SERVERS, p in PAIRS}: use_l[pair_left[p],s] >= map[p,s] ;
s.t. use_r_cells{s in SERVERS, (p,r) in prs}: use_r[r,s] >= map[p,s] ;

solve;

param balance{s in SERVERS} := sum{c in PAIRS} (map[c,s] * pair_weight[c]);
param comunica{s in SERVERS} := 
	sum{r in RCELLS} (use_r[r,s] * rcell_network[s,r]) +
	sum{l in LCELLS} (use_l[l,s] * lcell_network[s,l]);

for {c in PAIRS} {
	printf "%s ", c;
	for {s in SERVERS: map[c,s] > 0.9999999} {
	    printf " %d", s;
	}
	for {s in SERVERS: map[c,s] < 1 && map[c,s] > 0} {
	    printf " %d (%3.2f)", s, map[c,s];
	}
	printf "\n";
}
printf "\n";

printf "ser:";
for {p in PAIRS} {
	printf " %4s", p;
}
printf " :%11s %11s\n", "balance", "comum";

for {s in SERVERS} {
	printf "%3i: ", s;
	for{c in PAIRS} {
		printf "%3.2f ", map[c,s];
	}
	printf ": %10i %10i\n", balance[s], comunica[s];
}
printf "\n";

printf "\n+ Resultado: %i, max %i\n", z, x0;

param average := sum{s in SERVERS}(balance[s])/card(SERVERS);
param variance := sum{s in SERVERS}(balance[s]*balance[s])/card(SERVERS) - average*average;
param communication := sum{s in SERVERS} comunica[s];

printf "+ Communication: %i, Variance: %f\n\n", communication, variance;

end;
