__kernel void computeBodiesAcceleration(
    __global const float *in_qx,
    __global const float *in_qy,
    __global const float *in_qz,
    __global const float *in_m,
    __global float *out_ax,
    __global float *out_ay,
    __global float *out_az,
    const unsigned long n,
    const float soft_squared,
    const float g
)
{
    unsigned long i = get_global_id(0);
    if (i > n)
        return;

    for (unsigned long j = 0; j < n; j++) {
        const float rijx = in_qx[j] - in_qx[i];
        const float rijy = in_qy[j] - in_qy[i];
        const float rijz = in_qz[j] - in_qz[i];

        // compute the || rij ||² distance between body i and body j
        const float rij_squared = rijx * rijx + rijy * rijy + rijz * rijz;

        // compute the acceleration value between body i and body j: || ai || = G.mj / (|| rij ||² + e²)^{3/2}
        const float ai = g * in_m[j] / ((rij_squared + soft_squared) * sqrt(rij_squared + soft_squared)); // 5 flops

        // add the acceleration value into the acceleration vector: ai += || ai ||.rij
        out_ax[i] += ai * rijx;
        out_ay[i] += ai * rijy;
        out_az[i] += ai * rijz;

        barrier(CLK_GLOBAL_MEM_FENCE);
    }
}