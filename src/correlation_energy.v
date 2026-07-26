// ============================================================
// correlation_energy.v  --  Area-optimized (iterative MAC) FIXED
// ============================================================
module correlation_energy #(
    parameter K  = 8,
    parameter N  = 25,
    parameter Q  = 25,
    parameter AW = 8,
    parameter RW = 32,
    parameter EW = 64
)(
    input  wire                  clk,
    input  wire                  rst_n,
    input  wire                  start,
    input  wire [K*N*AW-1:0]     A_flat,
    input  wire [K*Q*RW-1:0]     Rr_flat,
    input  wire [K*Q*RW-1:0]     Ri_flat,
    input  wire [N-1:0]          mask,
    output reg  [N*EW-1:0]       energy_flat,
    output reg                   done
);
    localparam KK = K;
    localparam NN = N;
    localparam QQ = Q;

    reg busy;
    reg [31:0] nn, qq, kk;

    reg signed [63:0] cr, ci;
    reg signed [63:0] en_acc;

    wire signed [AW-1:0] a_val  = A_flat[(kk*N+nn)*AW +: AW];
    wire signed [RW-1:0] rr_val = Rr_flat[(kk*Q+qq)*RW +: RW];
    wire signed [RW-1:0] ri_val = Ri_flat[(kk*Q+qq)*RW +: RW];

    wire signed [63:0] next_cr = cr + $signed(a_val) * $signed(rr_val);
    wire signed [63:0] next_ci = ci + $signed(a_val) * $signed(ri_val);

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            done <= 0; energy_flat <= 0;
            busy <= 0; nn <= 0; qq <= 0; kk <= 0;
            cr <= 0; ci <= 0; en_acc <= 0;
        end else begin
            done <= 0;

            if (start && !busy) begin
                busy <= 1;
                nn   <= 0; qq <= 0; kk <= 0;
                cr <= 0; ci <= 0; en_acc <= 0;
            end else if (busy) begin
                if (mask[nn]) begin
                    energy_flat[nn*EW +: EW] <= 0;
                    nn <= nn + 1;
                    qq <= 0; kk <= 0;
                    cr <= 0; ci <= 0; en_acc <= 0;

                    if (nn == NN-1) begin
                        done <= 1; busy <= 0;
                    end
                end else begin
                    cr <= next_cr;
                    ci <= next_ci;

                    if (kk == KK-1) begin
                        en_acc <= en_acc + (next_cr*next_cr) + (next_ci*next_ci);
                        cr <= 0; ci <= 0;
                        kk <= 0;

                        if (qq == QQ-1) begin
                            energy_flat[nn*EW +: EW] <= en_acc + (next_cr*next_cr) + (next_ci*next_ci);
                            en_acc <= 0;
                            qq <= 0;
                            if (nn == NN-1) begin
                                done <= 1; busy <= 0;
                            end else nn <= nn + 1;
                        end else qq <= qq + 1;
                    end else kk <= kk + 1;
                end
            end
        end
    end
endmodule