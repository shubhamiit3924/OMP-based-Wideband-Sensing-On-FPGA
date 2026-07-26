// ============================================================
// mat_mul_AX.v  --  Area-optimized (iterative MAC) FIXED
// ============================================================
module mat_mul_AX #(
    parameter K  = 8,
    parameter N  = 25,
    parameter Q  = 25,
    parameter AW = 8,
    parameter XW = 8,
    parameter YW = 32
)(
    input  wire                      clk,
    input  wire                      rst_n,
    input  wire                      start,
    input  wire [K*N*AW-1:0]         A_flat,
    input  wire [N*Q*XW-1:0]         Xr_flat,
    input  wire [N*Q*XW-1:0]         Xi_flat,
    output reg  [K*Q*YW-1:0]         Yr_flat,
    output reg  [K*Q*YW-1:0]         Yi_flat,
    output reg                       done
);
    localparam KK = K;
    localparam NN = N;
    localparam QQ = Q;

    reg busy;
    reg [31:0] kk, qq, nn;

    reg signed [YW-1:0] acc_r, acc_i;

    wire signed [AW-1:0] a_val  = A_flat[(kk*NN+nn)*AW +: AW];
    wire signed [XW-1:0] xr_val = Xr_flat[(nn*QQ+qq)*XW +: XW];
    wire signed [XW-1:0] xi_val = Xi_flat[(nn*QQ+qq)*XW +: XW];

    wire signed [YW-1:0] next_acc_r = acc_r + (a_val * xr_val);
    wire signed [YW-1:0] next_acc_i = acc_i + (a_val * xi_val);

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            done <= 0; busy <= 0;
            kk <= 0; qq <= 0; nn <= 0;
            acc_r <= 0; acc_i <= 0;
            Yr_flat <= 0; Yi_flat <= 0;
        end else begin
            done <= 0;

            if (start && !busy) begin
                busy <= 1;
                kk <= 0; qq <= 0; nn <= 0;
                acc_r <= 0; acc_i <= 0;
            end else if (busy) begin
                acc_r <= next_acc_r;
                acc_i <= next_acc_i;

                if (nn == NN-1) begin
                    Yr_flat[(kk*QQ+qq)*YW +: YW] <= next_acc_r;
                    Yi_flat[(kk*QQ+qq)*YW +: YW] <= next_acc_i;

                    acc_r <= 0; acc_i <= 0;
                    nn <= 0;

                    if (qq == QQ-1) begin
                        qq <= 0;
                        if (kk == KK-1) begin
                            done <= 1;
                            busy <= 0;
                        end else kk <= kk + 1;
                    end else qq <= qq + 1;
                end else nn <= nn + 1;
            end
        end
    end
endmodule