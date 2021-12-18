+++
title = "A Strange but Effective Debugging Technique in Julia"
hascode = true
date = Date(2021, 12, 17)
rss = "This post introduces a strange but effective debugging technique in Julia involving global variables."
+++

@def tags = ["programming"]

Often, most of the work during debugging is to reconstruct out the intermediate
states in an algorithm. A debugger can often facilitate this need, but it could
be rather cumbersome to use if the breakpoint condition is complicated.

In Julia, we usually develop a package with a Julia [REPL](https://docs.julialang.org/en/v1/stdlib/REPL/)
open on the side. All the code executed in the REPL resides in the `Main`
module. Hence, one can access REPL variables via `Main.a_variable_name` in the
package one is developing. While debugging, we need to send the intermediate
states of a function to Julia REPL. Unfortunately, Julia forbids rebinding of a
variable from another module. One can bypass this limitation by defining an
unused variable like `_a = Ref{Any}()` in the REPL. Then, modify this variable
via `Main._a[] = a_state_of_interest_1, a_state_of_interest_2` in some function
of a package.

\section{An Example}
```julia
julia> function foo(n)
           a = b = 1
           for i in 1:n
               a, b = b, a + b
               if a <= 0
                   Main._a[] = a, b
               end
           end
           return a
       end
foo (generic function with 1 method)

julia> _a = Ref{Any}()
Base.RefValue{Any}(#undef)

julia> foo(100)
1298777728820984005

julia> _a
Base.RefValue{Any}((-2437933049959450366, 3736710778780434371))
```
