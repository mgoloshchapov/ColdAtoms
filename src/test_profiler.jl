using ProfileView, Profile

function profile_test_sort(n, len=100000)
    for i = 1:n
        list = []
        for _ in 1:len
            push!(list, rand())
        end
        sort!(list)
    end
end

@profview profile_test_sort(1)  # run once to trigger compilation (ignore this one)
@profview profile_test_sort(10) # measure runtime
ProfileView.view()