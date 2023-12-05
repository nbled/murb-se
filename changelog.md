cpu+naive: 20* n^2
(miguel) 467.5673 FPS, 2138.84 ms

cpu+optim : 
reduce the number of operation by computing each one once.
30 * (n^2 - n) / 2
result : 
(miguel) 588.4414 FPS, 1699.821 ms

changing layout: AoS => SoA
(nolan) 626.3245 FPS, 1596.6419999999998 ms
(miguel) 617.6135000000002 FPS, 1619.385 ms

remove redondant computation
(miguel) 864.9227000000001 FPS, 1156.262 ms

remove pow when doing a square
(miguel) 906.2089 FPS, 1105.008 ms

change loop condition
(miguel) 1035.8300000000004 FPS, 965.5302000000001 ms

switched from pow(a, 3/2) to a * sqrt(a)
(nolan) 2993.1580000000004 FPS, 334.124 ms
(miguel) 3134.459 FPS, 319.0527 ms

simd first version (naive)
(miguel) 6008.78 FPS, 166.423 ms (on x86_64)
