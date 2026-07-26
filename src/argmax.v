// ============================================================
// argmax.v  --  Area-optimized (iterative scan)
// ============================================================
module argmax #(
    parameter N  = 25,
    parameter EW = 64,
    parameter IW = 5
)(
    input  wire             clk,
    input  wire             rst_n,
    input  wire             start,
    input  wire [N*EW-1:0]  energy_flat,
    output reg  [IW-1:0]    max_idx,
    output reg  [EW-1:0]    max_val,
    output reg              done
);
    localparam NN = N;

    reg busy;
    reg [IW-1:0] idx;
    reg [EW-1:0] cur_max;
    reg [IW-1:0] cur_idx;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            max_idx <= 0; max_val <= 0; done <= 0;
            idx <= 0; cur_max <= 0; cur_idx <= 0; busy <= 0;
        end else begin
            done <= 0;

            if (start && !busy) begin
                busy    <= 1;
                idx     <= 0;
                cur_max <= energy_flat[0 +: EW];
                cur_idx <= 0;
            end else if (busy) begin
                if (energy_flat[idx*EW +: EW] > cur_max) begin
                    cur_max <= energy_flat[idx*EW +: EW];
                    cur_idx <= idx;
                end

                if (idx == NN-1) begin
                    max_idx <= cur_idx;
                    max_val <= cur_max;
                    done    <= 1;
                    busy    <= 0;
                end else begin
                    idx <= idx + 1;
                end
            end
        end
    end
endmodule