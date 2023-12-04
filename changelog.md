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