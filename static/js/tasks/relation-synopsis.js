
var document_pk = relation_task_settings.document_pk;
var doc_api_data = null;

$( document ).ready(function(document_pk) {

    $.getJSON('/task/relation/'+ relation_task_settings.document_pk +'/analysis/', function(document_api_data) {
        doc_api_data = document_api_data;


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

            this.style.background="#f7f7f9";
            this.style.cursor="pointer";
            var get_id = (this.id).split("-")[2];
            current_id = get_id;
            relation_pk = current_id;
            function show_results(document_pk, relation_pk, doc_api_data) {
                relation_specific_answers = _.filter(doc_api_data, {id: parseInt(relation_pk)})
                var answers = relation_specific_answers[0]['answers']
                var answer_counts = _.countBy( _.map(answers, function(x) { return x['answer']['id']; }) );
                var answer_text = {};
                var personal_ann = '';

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

                var chart = d3.select('#chart-' + relation_pk).style('width', '100%');
                var bar = chart.selectAll('div')
                  .data(data)
                  .enter()
                    .append('div')
                      .attr('class', 'bar-component')
                      .style('width', function(d) { return ((d['value']/max)*100) + '%'; } )
                      .style('background-color', function(d, i) { return color(i); });

                var bold = function(d, i) {
                  if( i == 1 && d['self'] ) {
                    return '<strong>';
                  } else if( i == 2 && d['self'] ) {
                    return '</strong>';
                  } else {
                    return '';
                  }
                }

                var list = d3.select('#chart-list-' + relation_pk);
                var bar = list.selectAll('li')
                  .data(data)
                  .enter()
                    .append('li')
                      .html(function(d, i) { return '<div class="box" style="background-color:'+color(i)+';"></div> <p>' + bold(d, 1) + ((d['value']/max)*100).toFixed() + '% – ' + d['label'] + bold(d, 2) + '</p>'; });
            }

            show_results(document_pk, relation_pk, doc_api_data);

            current_circle = this;
            $('#relation-' + previous_id).hide();
            $('#relation-' + current_id).show();
            click_flag = true;
        })
    });
});


