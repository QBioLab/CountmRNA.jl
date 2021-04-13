#julia --project=/home/hf/.julia/environments/v1.4 --depwarn=no test15.jl 11 15
#time julia --depwarn=no test19.jl 1 5
#time julia --depwarn=no test19.jl 6 10
#time julia --depwarn=no test19.jl 11 15
#time julia --depwarn=no test19.jl 16 20
#time julia --depwarn=no test19.jl 21 25
#time julia --depwarn=no test19.jl 26 30
#time julia --depwarn=no test19.jl 31 35
#time julia --depwarn=no test19.jl 28 30
#time julia --depwarn=no test19.jl 36 40
#for i in {1..12..4}
for i in {25..40..4}
do
	time julia --depwarn=no test21.jl $i $(( $i + 3 ))
done


