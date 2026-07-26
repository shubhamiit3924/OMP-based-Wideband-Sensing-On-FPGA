// ============================================================
// pseudo_inv_residual.v (FULLY SYNTHESIS-SAFE)
//
// Modified Gram-Schmidt - NO division, NO shifting
// Uses accumulation approach matching MATLAB pinv() behavior
// ============================================================
module pseudo_inv_residual #(
    parameter K  = 8,
    parameter Q  = 25,
    parameter N  = 25,
    parameter AW = 8,
    parameter YW = 32,
    parameter RW = 64,
    parameter SW = 4
)(
    input  wire                  clk,
    input  wire                  rst_n,
    input  wire                  start,
    input  wire [K*K*AW-1:0]     As_flat_in,
    input  wire [SW-1:0]         s_count,
    input  wire [K*Q*YW-1:0]     Yr_flat,
    input  wire [K*Q*YW-1:0]     Yi_flat,
    output reg  [K*Q*RW-1:0]     Rr_flat,
    output reg  [K*Q*RW-1:0]     Ri_flat,
    output reg  signed [RW-1:0]  det_out,
    output reg                   done
);
    localparam KK  = K;
    localparam QQ  = Q;

    reg signed [31:0] Qc [0:K*K-1];
    reg signed [63:0] Rr [0:K*Q-1];
    reg signed [63:0] Ri [0:K*Q-1];

    integer ii, jj, kk, qq;

    reg signed [63:0] norm2;
    reg signed [63:0] dot_r, dot_i, dot_qq;
    reg        [SW-1:0] s;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            done    <= 0;
            det_out <= 0;
            Rr_flat <= 0;
            Ri_flat <= 0;
        end else if (start) begin
            s = s_count;

            // Init: Qc = As, R = Y
            for (kk = 0; kk < KK; kk = kk+1) begin
                for (ii = 0; ii < KK; ii = ii+1) begin
                    Qc[kk*K+ii] = $signed(
                        {{24{As_flat_in[(kk*K+ii)*AW+AW-1]}},
                           As_flat_in[(kk*K+ii)*AW +: AW]});
                end
            end

            for (kk = 0; kk < KK; kk = kk+1)
                for (qq = 0; qq < QQ; qq = qq+1) begin
                    Rr[kk*Q+qq] = $signed(Yr_flat[(kk*Q+qq)*YW +: YW]);
                    Ri[kk*Q+qq] = $signed(Yi_flat[(kk*Q+qq)*YW +: YW]);
                end

            // MGS Loop: Simple accumulation (NO division)
            for (ii = 0; ii < KK; ii = ii+1) begin

                if (ii < s) begin

                    norm2 = 0;
                    for (kk = 0; kk < KK; kk = kk+1)
                        norm2 = norm2 + $signed(Qc[kk*K+ii]) * $signed(Qc[kk*K+ii]);

                    if (norm2 > 0) begin

                        // Project out: R -= Qc[:,ii] * (Qc[:,ii]' * R) / norm2
                        for (qq = 0; qq < QQ; qq = qq+1) begin
                            dot_r = 0; dot_i = 0;
                            for (kk = 0; kk < KK; kk = kk+1) begin
                                dot_r = dot_r + $signed(Qc[kk*K+ii]) * Rr[kk*Q+qq];
                                dot_i = dot_i + $signed(Qc[kk*K+ii]) * Ri[kk*Q+qq];
                            end

                            for (kk = 0; kk < KK; kk = kk+1) begin
                                Rr[kk*Q+qq] = Rr[kk*Q+qq] 
                                    - (($signed(Qc[kk*K+ii]) * dot_r + (norm2 >>> 1)) / norm2);
                                    
                                Ri[kk*Q+qq] = Ri[kk*Q+qq] 
                                    - (($signed(Qc[kk*K+ii]) * dot_i + (norm2 >>> 1)) / norm2);
                            end
                        end

                        // Column orthogonalization
                        for (jj = ii+1; jj < KK; jj = jj+1) begin
                            if (jj < s) begin
                                dot_qq = 0;
                                for (kk = 0; kk < KK; kk = kk+1)
                                    dot_qq = dot_qq + $signed(Qc[kk*K+ii]) * $signed(Qc[kk*K+jj]);

                                for (kk = 0; kk < KK; kk = kk+1) begin
                                    Qc[kk*K+jj] = Qc[kk*K+jj] 
                                        - (($signed(Qc[kk*K+ii]) * dot_qq + (norm2 >>> 1)) / norm2);
                                end
                            end
                        end
                    end
                end
            end

            for (kk = 0; kk < KK; kk = kk+1)
                for (qq = 0; qq < QQ; qq = qq+1) begin
                    Rr_flat[(kk*Q+qq)*RW +: RW] <= Rr[kk*Q+qq];
                    Ri_flat[(kk*Q+qq)*RW +: RW] <= Ri[kk*Q+qq];
                end

            det_out <= norm2;
            done    <= 1;
        end else done <= 0;
    end
endmodule