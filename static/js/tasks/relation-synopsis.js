


function show_results(document_pk, relation_pk) {
  $('#tree-action-area').hide();


  $.getJSON('/task/relation/'+ document_pk +'/analysis/', function(api_data) {
        console.log('data');
        console.log(api_data);
    });

  $.getJSON('/task/relation/'+ document_pk +'/analysis/' + relation_pk + '/', function(api_data) {
    var answers = api_data[0]['answers']

    /* obj, key = identifier, value = count */
    var answer_counts = _.countBy( _.map(answers, function(x) { return x['answer']['id']; }) );

    var answer_text = {};
    var personal_ann = '';
    console.log('lalalalala');

    _.each(answers, function(a) {
      answer_text[a.answer.id] = a.answer.text;
      if(a['self']) { personal_ann = a.answer.id; }
    });

    var data = [];
    _.each(_.keys(answer_counts), function(answer_key) {
      data.push({
        'id': answer_key,
        'value': answer_counts[answer_key],
        'label': answer_text[answer_key],
        'self': answer_key == personal_ann
      });
    });

    data = _.sortBy(data, function(d) { return -d['value'] });

    var max = d3.sum(data.map(function(i){ return i['value']; }));
    var color = d3.scale.category20c();

    var chart = d3.select('#chart').style('width', '100%');
    var bar = chart.selectAll('div')
      .data(data)
      .enter()
        .append('div')
          .attr('class', 'bar-component')
          .style('width', function(d) { return ((d['value']/max)*100) + '%'; } )
          .style('background-color', function(d, i) { return color(i); });

    var bold = function(d, i) {
      console.log(i, d);
      if( i == 1 && d['self'] ) {
        return '<strong>';
      } else if( i == 2 && d['self'] ) {
        return '</strong>';
      } else {
        return '';
      }
    }

    var list = d3.select('#chart-list');
    var bar = list.selectAll('li')
      .data(data)
      .enter()
        .append('li')
          .html(function(d, i) { return '<div class="box" style="background-color:'+color(i)+';"></div> <p>' + bold(d, 1) + ((d['value']/max)*100).toFixed() + '% – ' + d['label'] + bold(d, 2) + '</p>'; });

    $('#feedback-next-action-area').addClass('shown').slideDown('fast');
  });
}


// Click method option
$( ".relation-circle" ).hover(function () {
    this.style.cursor="pointer";
    });

var current_id = "";
var previous_id = "";
var previous_circle = "";
var current_circle = "";
var click_flag = false;

$( ".relation-circle" ).mousedown(function () {

    previous_id = current_id;
    previous_circle = current_circle;
    if (click_flag == true) {
        previous_circle.style.background="transparent";

    };

    $('#relation-' + previous_id).hide();
    this.style.background="#00ff00";
    this.style.cursor="pointer";
    var get_id = (this.id).split("-")[2];
    current_id = get_id;
    current_circle = this;
    console.log(current_circle);
    $('#relation-' + previous_id).hide();
    $('#relation-' + current_id).show();

    click_flag = true;

})
