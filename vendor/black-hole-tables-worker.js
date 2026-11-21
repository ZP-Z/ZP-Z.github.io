import { buildTables } from './black-hole-tables.js';

self.onmessage = ({ data: distance }) => {
  const tables = buildTables(distance);
  self.postMessage(tables, Object.values(tables).map(data => data.buffer));
};
