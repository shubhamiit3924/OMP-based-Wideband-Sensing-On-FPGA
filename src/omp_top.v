// ============================================================
// omp_top.v  --  Top-level OMP Module (area-optimized flow) FIXED
// ============================================================
module omp_top #(
    parameter K  = 8,
    parameter N  = 25,
    parameter Q  = 25,
    parameter AW = 8,
    parameter XW = 8,
    parameter YW = 32,
    parameter RW = 64,
    parameter EW = 64,
    parameter IW = 5,
    parameter SW = 4
)(
    input  wire                clk,
    input  wire                rst_n,
    input  wire                start,
    input  wire [K*N*AW-1:0]   A_flat,
    input  wire [N*Q*XW-1:0]   Xr_flat,
    input  wire [N*Q*XW-1:0]   Xi_flat,
    input  wire [EW-1:0]       epsilon_sq,
    output reg  [K*IW-1:0]     found_bands_flat,
    output reg  [SW-1:0]       num_bands,
    output reg                 done
);
    localparam KQ = K*Q;

    localparam [3:0]
        ST_IDLE        = 4'd0,
        ST_WAIT_MEAS   = 4'd1,
        ST_WAIT_FRO    = 4'd2,
        ST_WAIT_ENERGY = 4'd3,
        ST_WAIT_ARGMAX = 4'd4,
        ST_WAIT_PINV   = 4'd5,
        ST_DONE        = 4'd6;

    reg [3:0]    state;
    reg [SW-1:0] iter;

    reg  ax_start, fro_start, en_start, am_start, pi_start;
    wire ax_done, fro_done, en_done, am_done, pi_done;

    wire [K*Q*YW-1:0] Yr_flat, Yi_flat;
    wire [EW-1:0]     fro_sq_val;
    wire             fro_converged;
    wire [N*EW-1:0]   energy_flat;
    wire [IW-1:0]     max_idx;
    wire [EW-1:0]     max_val;

    // ===== Modules =====
    mat_mul_AX #(.K(K),.N(N),.Q(Q),.AW(AW),.XW(XW),.YW(YW)) u_ax (
        .clk(clk), .rst_n(rst_n), .start(ax_start),
        .A_flat(A_flat), .Xr_flat(Xr_flat), .Xi_flat(Xi_flat),
        .Yr_flat(Yr_flat), .Yi_flat(Yi_flat), .done(ax_done)
    );

    reg [K*Q*RW-1:0] Rr_flat_reg;
    reg [K*Q*RW-1:0] Ri_flat_reg;

    // --- Correct per-element truncation (RW -> YW) ---
    wire [K*Q*YW-1:0] Rr32_flat;
    wire [K*Q*YW-1:0] Ri32_flat;
    genvar gi;
    generate
        for (gi = 0; gi < K*Q; gi = gi+1) begin : g_trunc
            assign Rr32_flat[gi*YW +: YW] = Rr_flat_reg[gi*RW +: YW];
            assign Ri32_flat[gi*YW +: YW] = Ri_flat_reg[gi*RW +: YW];
        end
    endgenerate

    fro_norm_sq #(.K(K),.Q(Q),.RW(RW),.FSW(EW)) u_fro (
        .clk(clk), .rst_n(rst_n), .start(fro_start),
        .Rr_flat(Rr_flat_reg), .Ri_flat(Ri_flat_reg),
        .epsilon_sq(epsilon_sq),
        .fro_sq(fro_sq_val), .converged(fro_converged), .done(fro_done)
    );

    reg [N-1:0] band_mask;

    correlation_energy #(.K(K),.N(N),.Q(Q),.AW(AW),.RW(YW),.EW(EW)) u_energy (
        .clk(clk), .rst_n(rst_n), .start(en_start),
        .A_flat(A_flat),
        .Rr_flat(Rr32_flat),
        .Ri_flat(Ri32_flat),
        .mask(band_mask),
        .energy_flat(energy_flat), .done(en_done)
    );

    argmax #(.N(N),.EW(EW),.IW(IW)) u_argmax (
        .clk(clk), .rst_n(rst_n), .start(am_start),
        .energy_flat(energy_flat),
        .max_idx(max_idx), .max_val(max_val), .done(am_done)
    );

    wire [K*Q*RW-1:0] Rr_new_flat;
    wire [K*Q*RW-1:0] Ri_new_flat;
    wire signed [RW-1:0] det_out;

    reg [K*K*AW-1:0] As_flat_reg;
    reg [SW-1:0] s_count_reg;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) s_count_reg <= 0;
        else if (state == ST_WAIT_ARGMAX && am_done)
            s_count_reg <= iter + 1;
    end

    pseudo_inv_residual #(
        .K(K),.Q(Q),.N(N),.AW(AW),.YW(YW),.RW(RW),.SW(SW)
    ) u_pinv (
        .clk(clk), .rst_n(rst_n), .start(pi_start),
        .As_flat_in(As_flat_reg), .s_count(s_count_reg),
        .Yr_flat(Yr_flat), .Yi_flat(Yi_flat),
        .Rr_flat(Rr_new_flat), .Ri_flat(Ri_new_flat),
        .det_out(det_out), .done(pi_done)
    );

    integer ii;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state            <= ST_IDLE;
            iter             <= 0;
            num_bands        <= 0;
            band_mask        <= 0;
            done             <= 0;
            ax_start         <= 0;
            fro_start        <= 0;
            en_start         <= 0;
            am_start         <= 0;
            pi_start         <= 0;
            found_bands_flat <= 0;
            As_flat_reg      <= 0;
            Rr_flat_reg      <= 0;
            Ri_flat_reg      <= 0;
        end else begin
            ax_start <= 0; fro_start <= 0; en_start <= 0;
            am_start <= 0; pi_start  <= 0;

            case (state)
                ST_IDLE: begin
                    done <= 0;
                    if (start) begin
                        iter      <= 0;
                        num_bands <= 0;
                        band_mask <= 0;
                        ax_start  <= 1;
                        state     <= ST_WAIT_MEAS;
                    end
                end

                ST_WAIT_MEAS: begin
                    if (ax_done) begin
                        for (ii = 0; ii < K*Q; ii = ii+1) begin
                            Rr_flat_reg[ii*RW +: RW] <=
                                {{(RW-YW){Yr_flat[ii*YW+YW-1]}},
                                  Yr_flat[ii*YW +: YW]};
                            Ri_flat_reg[ii*RW +: RW] <=
                                {{(RW-YW){Yi_flat[ii*YW+YW-1]}},
                                  Yi_flat[ii*YW +: YW]};
                        end
                        fro_start <= 1;
                        state     <= ST_WAIT_FRO;
                    end
                end

                ST_WAIT_FRO: begin
                    if (fro_done) begin
                        if (fro_converged || iter >= K)
                            state <= ST_DONE;
                        else begin
                            en_start <= 1;
                            state    <= ST_WAIT_ENERGY;
                        end
                    end
                end

                ST_WAIT_ENERGY: begin
                    if (en_done) begin
                        am_start <= 1;
                        state    <= ST_WAIT_ARGMAX;
                    end
                end

                ST_WAIT_ARGMAX: begin
                    if (am_done) begin
                        found_bands_flat[iter*IW +: IW] <= max_idx;
                        band_mask[max_idx]              <= 1'b1;
                        for (ii = 0; ii < K; ii = ii+1)
                            As_flat_reg[(ii*K+iter)*AW +: AW] <=
                                A_flat[(ii*N+max_idx)*AW +: AW];
                        iter     <= iter + 1;
                        pi_start <= 1;
                        state    <= ST_WAIT_PINV;
                    end
                end

                ST_WAIT_PINV: begin
                    if (pi_done) begin
                        Rr_flat_reg <= Rr_new_flat;
                        Ri_flat_reg <= Ri_new_flat;
                        fro_start   <= 1;
                        state       <= ST_WAIT_FRO;
                    end
                end

                ST_DONE: begin
                    num_bands <= iter;
                    done      <= 1;
                    state     <= ST_IDLE;
                end

                default: state <= ST_IDLE;
            endcase
        end
    end
endmodule