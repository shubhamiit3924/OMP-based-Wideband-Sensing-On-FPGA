// ============================================================
// fro_norm_sq.v  --  Area-optimized (iterative accumulation) FIXED
// ============================================================
module fro_norm_sq #(
    parameter K   = 8,
    parameter Q   = 25,
    parameter RW  = 64,
    parameter FSW = 64
)(
    input  wire                   clk,
    input  wire                   rst_n,
    input  wire                   start,
    input  wire [K*Q*RW-1:0]      Rr_flat,
    input  wire [K*Q*RW-1:0]      Ri_flat,
    input  wire [FSW-1:0]         epsilon_sq,
    output reg  [FSW-1:0]         fro_sq,
    output reg                    converged,
    output reg                    done
);
    localparam KQ = K*Q;

    reg busy;
    reg [31:0] idx;
    reg [FSW-1:0] acc;

    wire signed [RW-1:0] rr = Rr_flat[idx*RW +: RW];
    wire signed [RW-1:0] ri = Ri_flat[idx*RW +: RW];

    wire [FSW-1:0] next_acc = acc
        + ($signed(rr) * $signed(rr))
        + ($signed(ri) * $signed(ri));

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            fro_sq <= 0; converged <= 0; done <= 0;
            acc <= 0; idx <= 0; busy <= 0;
        end else begin
            done <= 0;
            if (start && !busy) begin
                busy <= 1;
                idx  <= 0;
                acc  <= 0;
            end else if (busy) begin
                acc <= next_acc;

                if (idx == KQ-1) begin
                    fro_sq    <= next_acc;
                    converged <= (next_acc <= epsilon_sq);
                    done      <= 1;
                    busy      <= 0;
                end else idx <= idx + 1;
            end
        end
    end
endmodule