/*
 * $Id: flexi.h,v 1.2 2004/06/16 21:25:54 mbarnes Exp $
 */

/****************************************************************************
*  'FlexiLib' - Flexible Input Library                                      *
*  Copyright (C) 2000  Gerhard J. Theurich (gjt)                            *
*                                                                           *
*  This program is free software; you can redistribute it and/or modify     *
*  it under the terms of the GNU General Public License as published by     *
*  the Free Software Foundation; either version 2 of the License, or        *
*  (at your option) any later version.                                      *
*                                                                           *
*  This program is distributed in the hope that it will be useful,          *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
*  GNU General Public License for more details.                             *
*                                                                           *
*  You should have received a copy of the GNU General Public License        *
*  along with this program; if not, write to the Free Software              *
*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA  *
****************************************************************************/


int flexi_request_key(char *key, int required, char *fmt, ...);
  /* Request to add 'key' with format string 'fmt' into the current        */
  /* keyqueue. Returns 0 if no problem was detected and 1 if there were    */
  /* less arguments given than 'fmt' required.                             */


int flexi_gather_keys(FILE *fp);
  /* Gather all requested keys in current keyqueue. The corresponding      */
  /* arguments are set to those found at the last matching key in 'fp'.    */
  /* Returns the number of keys in current keyqueue that could not be      */
  /* resolved in 'fp'.                                                     */


int flexi_resolved_key(char *key);
  /* Returns 1 if 'key' was resolved by a previous 'flexi_gather_keys()'   */
  /* call, 0 if 'key' was not resolved and -1 if 'key' does not exist in   */
  /* current keyqueue.                                                     */
  
  
void flexi_clear_keys(void);
  /* Clear the entire current keyqueue.                                    */
