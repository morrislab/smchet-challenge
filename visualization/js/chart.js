google.load('visualization', '1.1', {packages: ['corechart']});

function get_json(url, callback) {
  var req = new XMLHttpRequest();
  req.open('GET', url, true);
  req.onreadystatechange = function (evt) {
    if (req.readyState == 4) {
       if(req.status == 200)
        callback(JSON.parse(req.responseText));
       else
        console.log('Error retrieving ' + url);
    }
  };
  req.send(null);
}

function prepare_frac_clonal_rows(ds1, ds2) {
  var rows = [];
  for(var sampid in ds1) {
    if(!(sampid in ds2)) {
      console.log('No data for ' + sampid + ' in ds2');
      continue;
    }

    if(ds1[sampid]['num_subclones'] > ds2[sampid]['num_subclones']) {
      var style = 'point { fill-color: #ffa500 }';
    } else if(ds1[sampid]['num_subclones'] < ds2[sampid]['num_subclones']) {
      var style = 'point { fill-color: #ff0000 }';
    } else {
      var style = '';
    }

    rows.push([
      ds1[sampid]['frac_clonal'],
      ds2[sampid]['frac_clonal'],
      sampid,
      style
    ]);
  }
  return rows;
}

function prepare_purity_rows(ds1, ds2) {
  var rows = [];
  for(var sampid in ds1) {
    if(!(sampid in ds2)) {
      console.log('No data for ' + sampid + ' in ds2');
      continue;
    }
    rows.push([ds1[sampid]['purity'], ds2[sampid]['purity'], sampid]);
  }
  return rows;
}

function prepare_subclone_counts_rows(ds1, ds2) {
  var aggregated = {};
  var dataset_names = {};
  var total_datasets = 0;

  for(var sampid in ds1) {
    if(!(sampid in ds2)) {
      console.log('No data for ' + sampid + ' in ds2');
      continue;
    }

    var key = ds1[sampid]['num_subclones'] + ',' + ds2[sampid]['num_subclones'];
    if(typeof aggregated[key] === 'undefined')
      aggregated[key] = 0;
    if(typeof dataset_names[key] === 'undefined')
      dataset_names[key] = [];

    total_datasets++;
    aggregated[key]++;
    dataset_names[key].push(sampid);
  }

  var max_size = 80;
  var rows = [];
  for(var key in aggregated) {
    var coords = key.split(',');
    var agg_count = aggregated[key];
    var pct = round(100 * agg_count / total_datasets, 1) + '%';
    dataset_names[key].sort();
    rows.push([
      parseInt(coords[0], 10),
      parseInt(coords[1], 10),
      'Number of datasets: ' + agg_count + ' (' + pct + ')\n\n' + dataset_names[key].join('\n'),
      'point { size: ' + (max_size * agg_count / total_datasets) + '}'
    ]);
  }
  return rows;
}

function build_ticks(vals) {
  var min = Math.floor(Math.min.apply(Math, vals));
  var max = Math.ceil(Math.max.apply(Math, vals));
  var ticks = [];
  for(var i = min; i <= max; i++) {
    ticks.push(i);
  }
  return ticks;
}

function draw_charts(ds1, ds2) {
  var rows = {
    purity: prepare_purity_rows(ds1, ds2),
    frac_clonal: prepare_frac_clonal_rows(ds1, ds2),
    num_subclones: prepare_subclone_counts_rows(ds1, ds2),
  };
  draw_chart('Cellularity (' + rows.purity.length + ' datasets)', 'purity', rows.purity);
  draw_chart('Proportion of clonal SSMs (' + rows.frac_clonal.length + ' datasets)', 'frac_clonal', rows.frac_clonal, {
    extra_cols: [{type: 'string', role: 'style'}]
  });

  var hvals = [], vvals = [];
  rows.num_subclones.forEach(function(row) {
    hvals.push(row[0]);
    vvals.push(row[1]);
  });
  draw_chart('Number of subclones', 'num_subclones', rows.num_subclones, {
    hticks: build_ticks(hvals),
    vticks: build_ticks(vvals),
    extra_cols: [{type: 'string', role: 'style'}]
  });

  draw_num_subclones_table(ds1, ds2);
}

function draw_num_subclones_table(ds1, ds2) {
  var table = document.getElementById('num_subclones_table');
  table.innerHTML += '<tr><th>Dataset</th><th>' + get_url_param('ds1n') + '</th><th>' + get_url_param('ds2n') + '</th></tr>';


  var rows = {
    gt: [],
    lt: [],
    equal: []
  };
  var bgcolor = {
    gt: '#ffd484',   // Orange
    lt: '#ff8484',   // Red
    equal: '#84d5ff' // Blue
  };

  for(var sampid in ds1) {
    if(!(sampid in ds2)) {
      console.log('No data for ' + sampid + ' in ds2');
      continue;
    }

    var cells = [sampid, ds1[sampid]['num_subclones'], ds2[sampid]['num_subclones']];
    if(ds1[sampid]['num_subclones'] > ds2[sampid]['num_subclones']) {
      rows.gt.push(cells);
    } else if(ds1[sampid]['num_subclones'] < ds2[sampid]['num_subclones']) {
      rows.lt.push(cells);
    } else {
      rows.equal.push(cells);
    }
  }

  ['equal', 'gt', 'lt'].forEach(function(relation) {
    rows[relation].sort(function(a, b) {
      return a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0;
    });
    rows[relation] = rows[relation].map(function(c) {
      var str = '<tr style="background-color: ' + bgcolor[relation] + '">';
      str += c.map(function(d) {
        return '<td>' + d + '</td>';
      }).join('');
      str += '</tr>';
      return str;
    });
    table.innerHTML += rows[relation].join('\n');
  });
}

function draw_chart(title, container_id, rows, options) {
  if(typeof options === 'undefined')
    options = {};

  var data = new google.visualization.DataTable();
  data.addColumn('number', 'ds1');
  data.addColumn('number', 'ds2');
  data.addColumn({type: 'string', role: 'tooltip'});

  if(options.extra_cols) {
    options.extra_cols.forEach(function(extra_col) {
      data.addColumn(extra_col);
    });
  }
  data.addRows(rows);

  var plot_options = {
    title: title,
    legend: 'none',
    width: 900,
    height: 500,
    hAxis: { title: get_url_param('ds1n') },
    vAxis: { title: get_url_param('ds2n') },
  };
  if(options.hticks) {
    plot_options.hAxis.ticks = options.hticks;
  }
  if(options.vticks) {
    plot_options.vAxis.ticks = options.vticks;
  }

  var chart = new google.visualization.ScatterChart(document.getElementById(container_id));
  chart.draw(data, plot_options);
}

function round(value, exp) {
  if (typeof exp === 'undefined' || +exp === 0)
    return Math.round(value);

  value = +value;
  exp  = +exp;

  if (isNaN(value) || !(typeof exp === 'number' && exp % 1 === 0))
    return NaN;

  // Shift
  value = value.toString().split('e');
  value = Math.round(+(value[0] + 'e' + (value[1] ? (+value[1] + exp) : exp)));

  // Shift back
  value = value.toString().split('e');
  return +(value[0] + 'e' + (value[1] ? (+value[1] - exp) : -exp));
}

function get_url_param(name) {
  name = name.replace(/[\[]/, "\\[").replace(/[\]]/, "\\]");
  var regex = new RegExp("[\\?&]" + name + "=([^&#]*)"),
      results = regex.exec(location.search);
  return results === null ? "" : decodeURIComponent(results[1].replace(/\+/g, " "));
}

function main() {
  get_json('data/' + get_url_param('ds1') + '.json', function(ds1) {
    get_json('data/' + get_url_param('ds2') + '.json', function(ds2) {
      document.title = get_url_param('title');
      draw_charts(ds1, ds2);
    });
  });
}

google.setOnLoadCallback(main);
