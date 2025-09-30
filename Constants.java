/*
 * Copyright (c) 2015, 2019, 2020, 2023, 2025 Russell A. Brown
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/**
 * @author Russell A. Brown
 */
	
/**
 * <p>
 * The Constants class specifies various constants.
 * <p>
 */
public class Constants {

	protected static final int INSERTION_SORT_CUTOFF = 15;
		
    protected static final int MEDIAN_OF_MEDIANS_CUTOFF = 15;

    protected static final int NLOGN_CUTOFF = 4096;

    protected static final int KNLOGN_CUTOFF = 4096;

    protected static final int MERGE_CUTOFF = 4096;

    protected static final long MILLISECONDS_TO_SECONDS = 1000L;

    protected static final int MULTI_THREAD_CUTOFF = 65536;

    protected static final int HEIGHT_DIFF = 1;
    
    public static final boolean AVL_BALANCE = false;

    public static final boolean NLOGN = true;

    public static final boolean KD_MAP_DYNAMIC = true;

    public static final boolean ENABLE_1TO3 = true;

    public static final boolean ENABLE_PREFERRED_NODE = true;

    public static final boolean ENABLE_LINKED_LIST = true;

} // class Constants
