// ============================================================
// omp_top_bram.v  --  Memory-mapped BRAM interface for OMP (FIXED)
// ============================================================
// Replaces flat I/O with addressable memory interface
// ~20 I/O pins instead of 29,000+ bits
// ============================================================

module omp_top_bram #(
    parameter K  = 8,
    parameter N  = 25,
    parameter Q  = 25,
    parameter AW = 8,
    parameter XW = 8,
    parameter YW = 32,
    parameter RW = 64,
    parameter EW = 64,
    parameter IW = 5,
    parameter SW = 4,
    parameter ADDR_WIDTH = 16,
    parameter DATA_WIDTH = 32
)(
    input  wire                      clk,
    input  wire                      rst_n,
    
    // Memory-mapped interface (~20 pins)
    input  wire [ADDR_WIDTH-1:0]     mem_addr,
    input  wire [DATA_WIDTH-1:0]     mem_data_in,
    output wire [DATA_WIDTH-1:0]     mem_data_out,
    input  wire                      mem_write_en,
    input  wire                      mem_read_en,
    
    // Status outputs
    output wire [SW-1:0]             num_bands_out,
    output wire                      done_out
);

    // ===== Memory Layout (Simplified) =====
    // 0x0000-0x07FF: A_matrix (input)
    // 0x0800-0x1FFF: Xr_matrix (input)
    // 0x2000-0x39FF: Xi_matrix (input)
    // 0x3A00-0x51FF: Yr_matrix (output)
    // 0x5200-0x69FF: Yi_matrix (output)
    // 0x6A00-0x6A27: found_bands (output, 1 index per 32-bit word)
    // 0x6A28-0x6A2B: num_bands (output)
    // 0x6A2C-0x6A2F: control/status register

    localparam A_START     = 16'h0000;
    localparam A_END       = 16'h07FF;
    localparam XR_START    = 16'h0800;
    localparam XR_END      = 16'h1FFF;
    localparam XI_START    = 16'h2000;
    localparam XI_END      = 16'h39FF;
    localparam YR_START    = 16'h3A00;
    localparam YR_END      = 16'h51FF;
    localparam YI_START    = 16'h5200;
    localparam YI_END      = 16'h69FF;
    localparam FB_START    = 16'h6A00;
    localparam FB_END      = 16'h6A27;
    localparam NB_ADDR     = 16'h6A28;
    localparam CTRL_ADDR   = 16'h6A2C;

    // ===== Internal BRAM Arrays (32-bit words) =====
    reg [DATA_WIDTH-1:0] bram_A   [0:49];     // A_matrix (50 × 32 bits)
    reg [DATA_WIDTH-1:0] bram_Xr  [0:156];    // Xr_matrix (157 × 32 bits)
    reg [DATA_WIDTH-1:0] bram_Xi  [0:156];    // Xi_matrix (157 × 32 bits)
    reg [DATA_WIDTH-1:0] bram_Yr  [0:199];    // Yr_matrix (200 × 32 bits) [OUTPUT]
    reg [DATA_WIDTH-1:0] bram_Yi  [0:199];    // Yi_matrix (200 × 32 bits) [OUTPUT]
    reg [DATA_WIDTH-1:0] bram_fb  [0:7];      // found_bands (8 × 32 bits = 40 bits) [OUTPUT]
    reg [DATA_WIDTH-1:0] num_bands_reg;
    reg [DATA_WIDTH-1:0] ctrl_reg;

    // ===== Intermediate flat vectors (reg type for procedural assignment) =====
    reg [K*N*AW-1:0]   A_flat;
    reg [N*Q*XW-1:0]   Xr_flat;
    reg [N*Q*XW-1:0]   Xi_flat;
    reg [K*Q*YW-1:0]   Yr_flat;
    reg [K*Q*YW-1:0]   Yi_flat;
    reg [K*IW-1:0]     found_bands_flat;

    // ===== Convert BRAM to flat vectors (combinatorial) =====
    genvar gi;
    generate
        // A_matrix: read 1600 bits from 50 × 32-bit words
        for (gi = 0; gi < K*N; gi = gi+1) begin : g_a_conv
            always @(*) begin
                A_flat[gi*AW +: AW] = 
                    bram_A[gi*AW/32][(gi*AW) % 32 +: AW];
            end
        end
        
        // Xr_matrix: read 5000 bits from 157 × 32-bit words
        for (gi = 0; gi < N*Q; gi = gi+1) begin : g_xr_conv
            always @(*) begin
                Xr_flat[gi*XW +: XW] = 
                    bram_Xr[gi*XW/32][(gi*XW) % 32 +: XW];
            end
        end
        
        // Xi_matrix: read 5000 bits from 157 × 32-bit words
        for (gi = 0; gi < N*Q; gi = gi+1) begin : g_xi_conv
            always @(*) begin
                Xi_flat[gi*XW +: XW] = 
                    bram_Xi[gi*XW/32][(gi*XW) % 32 +: XW];
            end
        end
    endgenerate

    // ===== Original OMP Core (unchanged) =====
    wire ax_done, fro_done, en_done, am_done, pi_done;
    reg  ax_start, fro_start, en_start, am_start, pi_start;
    
    wire [N*EW-1:0]   energy_flat;
    wire [IW-1:0]     max_idx;
    wire [EW-1:0]     max_val;
    wire [EW-1:0]     fro_sq_val;
    wire              fro_converged;
    wire [K*Q*RW-1:0] Rr_flat_new, Ri_flat_new;
    wire signed [RW-1:0] det_out;

    // Output wires from computation modules
    wire [K*Q*YW-1:0] Yr_flat_comp;
    wire [K*Q*YW-1:0] Yi_flat_comp;

    reg [K*Q*RW-1:0] Rr_flat_reg, Ri_flat_reg;
    reg [N-1:0] band_mask;
    reg [K*K*AW-1:0] As_flat_reg;
    reg [SW-1:0] s_count_reg;

    wire [K*Q*YW-1:0] Rr32_flat, Ri32_flat;
    generate
        for (gi = 0; gi < K*Q; gi = gi+1) begin : g_rtrunc
            assign Rr32_flat[gi*YW +: YW] = Rr_flat_reg[gi*RW +: YW];
            assign Ri32_flat[gi*YW +: YW] = Ri_flat_reg[gi*RW +: YW];
        end
    endgenerate

    // Instantiate original computation modules
    mat_mul_AX #(.K(K),.N(N),.Q(Q),.AW(AW),.XW(XW),.YW(YW)) u_ax (
        .clk(clk), .rst_n(rst_n), .start(ax_start),
        .A_flat(A_flat), .Xr_flat(Xr_flat), .Xi_flat(Xi_flat),
        .Yr_flat(Yr_flat_comp), .Yi_flat(Yi_flat_comp), .done(ax_done)
    );

    fro_norm_sq #(.K(K),.Q(Q),.RW(RW),.FSW(EW)) u_fro (
        .clk(clk), .rst_n(rst_n), .start(fro_start),
        .Rr_flat(Rr_flat_reg), .Ri_flat(Ri_flat_reg),
        .epsilon_sq(64'd361),
        .fro_sq(fro_sq_val), .converged(fro_converged), .done(fro_done)
    );

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

    pseudo_inv_residual #(.K(K),.Q(Q),.N(N),.AW(AW),.YW(YW),.RW(RW),.SW(SW)) u_pinv (
        .clk(clk), .rst_n(rst_n), .start(pi_start),
        .As_flat_in(As_flat_reg), .s_count(s_count_reg),
        .Yr_flat(Yr_flat), .Yi_flat(Yi_flat),
        .Rr_flat(Rr_flat_new), .Ri_flat(Ri_flat_new),
        .det_out(det_out), .done(pi_done)
    );

    // ===== State Machine (unchanged from original) =====
    localparam [3:0]
        ST_IDLE        = 4'd0,
        ST_WAIT_MEAS   = 4'd1,
        ST_WAIT_FRO    = 4'd2,
        ST_WAIT_ENERGY = 4'd3,
        ST_WAIT_ARGMAX = 4'd4,
        ST_WAIT_PINV   = 4'd5,
        ST_DONE        = 4'd6;

    reg [3:0] state;
    reg [SW-1:0] iter;
    reg process_start;

    integer ii, jj;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state            <= ST_IDLE;
            iter             <= 0;
            process_start    <= 0;
            num_bands_reg    <= 0;
            band_mask        <= 0;
            ax_start         <= 0;
            fro_start        <= 0;
            en_start         <= 0;
            am_start         <= 0;
            pi_start         <= 0;
            As_flat_reg      <= 0;
            Rr_flat_reg      <= 0;
            Ri_flat_reg      <= 0;
            ctrl_reg         <= 0;
            Yr_flat          <= 0;
            Yi_flat          <= 0;
            found_bands_flat <= 0;
        end else begin
            ax_start <= 0; fro_start <= 0; en_start <= 0;
            am_start <= 0; pi_start  <= 0;

            case (state)
                ST_IDLE: begin
                    if (process_start) begin
                        iter      <= 0;
                        num_bands_reg <= 0;
                        band_mask <= 0;
                        ax_start  <= 1;
                        state     <= ST_WAIT_MEAS;
                        process_start <= 0;
                    end
                end

                ST_WAIT_MEAS: begin
                    if (ax_done) begin
                        // Capture output from mat_mul_AX
                        Yr_flat <= Yr_flat_comp;
                        Yi_flat <= Yi_flat_comp;
                        
                        for (ii = 0; ii < K*Q; ii = ii+1) begin
                            Rr_flat_reg[ii*RW +: RW] <=
                                {{(RW-YW){Yr_flat_comp[ii*YW+YW-1]}},
                                  Yr_flat_comp[ii*YW +: YW]};
                            Ri_flat_reg[ii*RW +: RW] <=
                                {{(RW-YW){Yi_flat_comp[ii*YW+YW-1]}},
                                  Yi_flat_comp[ii*YW +: YW]};
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
                        s_count_reg <= iter + 1;
                        pi_start <= 1;
                        state    <= ST_WAIT_PINV;
                    end
                end

                ST_WAIT_PINV: begin
                    if (pi_done) begin
                        Rr_flat_reg <= Rr_flat_new;
                        Ri_flat_reg <= Ri_flat_new;
                        fro_start   <= 1;
                        state       <= ST_WAIT_FRO;
                    end
                end

                ST_DONE: begin
                    num_bands_reg <= iter;
                    state         <= ST_IDLE;
                end

                default: state <= ST_IDLE;
            endcase
        end
    end

    // ===== Memory Read/Write Controller =====
    integer bram_idx;
    
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            // Initialize BRAM to zero
            for (bram_idx = 0; bram_idx < 200; bram_idx = bram_idx+1) begin
                bram_A[bram_idx]   <= 0;
                bram_Xr[bram_idx]  <= 0;
                bram_Xi[bram_idx]  <= 0;
                bram_Yr[bram_idx]  <= 0;
                bram_Yi[bram_idx]  <= 0;
            end
            for (bram_idx = 0; bram_idx < 8; bram_idx = bram_idx+1)
                bram_fb[bram_idx]  <= 0;
            num_bands_reg <= 0;
            ctrl_reg      <= 0;
        end else begin
            // WRITE Operations
            if (mem_write_en) begin
                if (mem_addr >= A_START && mem_addr <= A_END) begin
                    bram_A[mem_addr[7:2]] <= mem_data_in;
                end else if (mem_addr >= XR_START && mem_addr <= XR_END) begin
                    bram_Xr[(mem_addr - XR_START) >> 2] <= mem_data_in;
                end else if (mem_addr >= XI_START && mem_addr <= XI_END) begin
                    bram_Xi[(mem_addr - XI_START) >> 2] <= mem_data_in;
                end else if (mem_addr == CTRL_ADDR) begin
                    ctrl_reg <= mem_data_in;
                    if (mem_data_in[0]) process_start <= 1;  // Start signal
                end
            end
            
            // Update output BRAM from internal signals (every cycle)
            for (bram_idx = 0; bram_idx < 200; bram_idx = bram_idx+1) begin
                if (bram_idx < 200) begin
                    bram_Yr[bram_idx] <= Yr_flat[bram_idx*DATA_WIDTH +: DATA_WIDTH];
                    bram_Yi[bram_idx] <= Yi_flat[bram_idx*DATA_WIDTH +: DATA_WIDTH];
                end
            end
            
            // Update found_bands BRAM (one index per 32-bit word, zero-padded)
            for (bram_idx = 0; bram_idx < K; bram_idx = bram_idx+1) begin
                bram_fb[bram_idx] <= {{(DATA_WIDTH - IW){1'b0}}, found_bands_flat[bram_idx*IW +: IW]};
            end
        end
    end

    // ===== Memory Read Multiplexer =====
    reg [DATA_WIDTH-1:0] read_data;
    
    always @(*) begin
        read_data = 32'd0;
        if (mem_read_en) begin
            if (mem_addr >= A_START && mem_addr <= A_END) begin
                read_data = bram_A[mem_addr[7:2]];
            end else if (mem_addr >= XR_START && mem_addr <= XR_END) begin
                read_data = bram_Xr[(mem_addr - XR_START) >> 2];
            end else if (mem_addr >= XI_START && mem_addr <= XI_END) begin
                read_data = bram_Xi[(mem_addr - XI_START) >> 2];
            end else if (mem_addr >= YR_START && mem_addr <= YR_END) begin
                read_data = bram_Yr[(mem_addr - YR_START) >> 2];
            end else if (mem_addr >= YI_START && mem_addr <= YI_END) begin
                read_data = bram_Yi[(mem_addr - YI_START) >> 2];
            end else if (mem_addr >= FB_START && mem_addr <= FB_END) begin
                // Return one band index per 32-bit read
                read_data = bram_fb[(mem_addr - FB_START) >> 2];
            end else if (mem_addr == NB_ADDR) begin
                read_data = num_bands_reg;
            end else if (mem_addr == CTRL_ADDR) begin
                read_data = {31'd0, (state == ST_DONE)};  // done flag in LSB
            end
        end
    end

    assign mem_data_out = read_data;
    assign num_bands_out = num_bands_reg[SW-1:0];
    assign done_out = (state == ST_DONE);

endmodule